import os
from src.cvsegmenter import CVSegmenter
from src.cvstitch import CVMaskStitcher
from src.cvmask import CVMask
from src import cvutils
from src import cvvisualize
from src import fcswrite
import src.walkthrough_utils as wu
from PIL import Image
import skimage
import numpy as np
import pandas as pd
import tensorflow as tf
from collections import defaultdict
import time
from matplotlib.pyplot import imshow, show
from tifffile import imsave

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

INPUT_PATH = "/mnt/codex/B005_CL"
OUTPUT_PATH_NAME = "/mnt/codex/B005_CL_segmentation"


VALID_IMAGE_EXTENSIONS = ('tif', 'jpg', 'png')

###########################



FILENAMES = [f for f in os.listdir(INPUT_PATH) if f.endswith(
            VALID_IMAGE_EXTENSIONS) and not f.startswith('.') ]
if len(FILENAMES) < 1:
    raise NameError(
        'No image files found.  Make sure you are pointing to the right directory')
else:
    print("Available images to segment at " + os.path.join(INPUT_PATH, 'bestFocus'))
    print(FILENAMES)


IS_CODEX_OUTPUT = True

NUCLEAR_CHANNEL_NAME = "Hoechst1"

###########################



CHANNEL_NAMES = []
ext = FILENAMES[0].split('.')[-1]
if 'tif' in ext:
    CHANNEL_NAMES = pd.read_csv(
        os.path.join(INPUT_PATH, "channelNames.txt"), sep='\t', header=None).values[:, 0]
    if NUCLEAR_CHANNEL_NAME not in CHANNEL_NAMES:
        raise NameError(
            'Nuclear channel name not found. Please double check channelNames.txt file')
    print('Picking channel', NUCLEAR_CHANNEL_NAME, 'from',
        len(CHANNEL_NAMES), 'total to segment on')
    print('Channel names:')
    print(CHANNEL_NAMES)


reference_image_path = os.path.join(INPUT_PATH, FILENAMES[0])

N_DIMS, EXT, DTYPE, SHAPE, READ_METHOD = cvutils.meta_from_image(reference_image_path)

print("Working with images of shape:", SHAPE)

OVERLAP = 80
THRESHOLD = 20
INCREASE_FACTOR = 2.5

stitcher = CVMaskStitcher(overlap=OVERLAP)

segmenter = CVSegmenter(
        SHAPE,
        os.path.join("src", "modelFiles", "final_weights.h5"),
        OVERLAP,
        INCREASE_FACTOR,
        THRESHOLD
    )
rows, cols = None, None
dataframe_regs = defaultdict(list)
columns = []
path = ''

print("\n", "Names of images to segment")
print(FILENAMES)

for CURR_IM_NAME in FILENAMES[3:]:

	path = os.path.join(INPUT_PATH, CURR_IM_NAME)
	image = np.array(READ_METHOD(path))

	if 'tif' in ext:
	    if N_DIMS == 4:
	        image = np.transpose(image, (2, 3, 0, 1))
	    elif N_DIMS == 3:
	        image = np.transpose(image, (1, 2, 0))
	image = image.reshape(SHAPE)
	nuclear_index = None
	if 'tif' in ext:
	    nuclear_index = cvutils.get_channel_index(NUCLEAR_CHANNEL_NAME, CHANNEL_NAMES)
	nuclear_image = cvutils.get_nuclear_image(N_DIMS-1, image, nuclear_index=nuclear_index)

	print("Segmenting " + CURR_IM_NAME + " on channel " + NUCLEAR_CHANNEL_NAME)
	masks, rows, cols = segmenter.segment_image(nuclear_image)
	#stitched_mask = stitcher.stitch_masks(masks, rows, cols)
	print("Stitching together cropped segmented subimages")
	stitched_mask = CVMask(stitcher.stitch_masks(masks, rows, cols))
	instances = stitched_mask.n_instances()
	print(instances, 'cell masks found by segmenter')
	if instances == 0:
	    print('No cells found in', filename)

	GROW_MASKS = True

	GROWTH_PIXELS = 3

	GROWTH_METHOD = 'Sequential'

	###########################



	if GROW_MASKS:
	    print("Computing centroids and bounding boxes for the masks.")
	    stitched_mask.compute_centroids()
	    stitched_mask.compute_boundbox()
	    print(f"Growing masks by {GROWTH_PIXELS} pixels")
	    stitched_mask.grow_masks(GROWTH_PIXELS, GROWTH_METHOD)
	
	SHOULD_COMPENSATE = True

	CSV_OR_FCS = 'csv'

	###########################

	stitched_mask.compute_centroids()
	stitched_mask.compute_boundbox()

	output_masks = stitched_mask

	QUANTIFICATION_OUTPUT_PATH = os.path.join(OUTPUT_PATH_NAME, 'quantifications')

	if not os.path.isdir(QUANTIFICATION_OUTPUT_PATH):
	    os.makedirs(QUANTIFICATION_OUTPUT_PATH)

	wu.compute_stats(output_masks, CURR_IM_NAME, image, IS_CODEX_OUTPUT, CHANNEL_NAMES, 
	                 SHOULD_COMPENSATE, QUANTIFICATION_OUTPUT_PATH, CSV_OR_FCS, GROWTH_PIXELS)


	VISUAL_OUTPUT_PATH = os.path.join(OUTPUT_PATH_NAME, 'visual_output')
	try:
	    os.makedirs(VISUAL_OUTPUT_PATH)
	except FileExistsError:
	    print("Directory already exists")

	print('Creating visual overlay output saved to', VISUAL_OUTPUT_PATH)

	new_path = os.path.join(VISUAL_OUTPUT_PATH, CURR_IM_NAME[:-4]) + '_mask.tif'

	outlines = cvvisualize.generate_mask_outlines(stitched_mask.flatmasks)
	imsave(new_path, outlines)
