
slib:
	cd polyeval/btfy/gf264 && $(MAKE) slib
	mv polyeval/btfy/gf264/libgf264btfy.so libgf264fft/

clean:
	cd polyeval/btfy/gf264 && $(MAKE) clean
	-rm libgf264fft/libgf264btfy.so
