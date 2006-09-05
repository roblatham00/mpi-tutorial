ALL: all-redirect

MAKE = make

all-redirect:
	cd examples && $(MAKE)

notes:
	cd examples && $(MAKE) notes

clean:
	cd examples && $(MAKE) clean
