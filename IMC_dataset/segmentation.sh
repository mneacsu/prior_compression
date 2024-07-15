cd steinbock

alias steinbock="docker run -v ./:/data -u $(id -u):$(id -g) --network host -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY ghcr.io/bodenmillergroup/steinbock"
steinbock segment deepcell --minmax
steinbock measure intensities
steinbock measure regionprops
# steinbock measure neighbors --type centroids --dmax 15

cd ..