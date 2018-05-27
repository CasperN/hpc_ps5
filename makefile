
info :
	@ echo "Usage:"
	@ echo "To see results: 'make go'"
	@ echo "To make raytracer: 'make rt'"
	@ echo "To run ray tracer: './rt 5000000 1000'"
	@ echo ""
	@ echo "For GPU version: 'make grt'"


go: rt
	./rt 5000000 1000
	python plotter.py result.bin

Readme.pdf : readme.md
	pandoc readme.md -o Readme.pdf

clean :
	rm rt result.bin

rt : ray_trace.c
	gcc -O3 -lm -Wall -std=c99 ray_trace.c -o rt

grt : gray_trace.cu
	nvcc -O3 -lm -Wall -std=c99 gray_trace.cu -o grt
