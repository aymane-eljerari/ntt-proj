docker run -itd \
	--name ntt-proj \
	--gpus all \
	-v $(pwd):/opt/ntt-proj \
	--privileged \
	aej-ntt-cuda

# disable rocm-gpu discovery
# --device=/dev/kfd \
# --device=/dev/dri \
# --group-add=video \
