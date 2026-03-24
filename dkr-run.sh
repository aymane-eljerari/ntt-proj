docker run -itd \
	--name ntt-proj \
	--device=/dev/kfd \
	--device=/dev/dri \
	--group-add=video \
	-v $(pwd):/opt/ntt-proj \
	aej-ntt-proj 
