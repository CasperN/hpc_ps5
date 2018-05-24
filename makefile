# CC = gcc -std=c99 -cc=/opt/local/bin/gcc-mp-7 -O3 -lm -Wall
CC = gcc -O3 -lm -Wall -std=c99

info :
	@ echo "Usage:"
	@ echo "To see results: 'make go'"
	@ echo "To make raytracer: 'make rt'"
	@ echo "To run ray tracer: './rt 1000 10000'"


go: rt
	./rt 5000000 1000
	python plotter.py result.bin

Readme.pdf : readme.md
	pandoc readme.md -o Readme.pdf

clean :
	rm rt result.bin

rt : main.c
	$(CC) main.c -o rt
