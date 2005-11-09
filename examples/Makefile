CTOPS=../maint/ctops

DIRS = life sparsemat pneo mutex 

PSFILES=life/mlife.ps \
	life/mlife-pt2pt.ps \
	life/mlife-io-stdout.ps \
	life/mlife-io-mpiio.ps \
	life/mlife-io-pnetcdf.ps \
	life/mlife-fence.ps \
	life/mlife-pscw.ps \
	life/mlife2d.ps \
	life/mlife2d-pt2pt.ps \
	life/mlife2d-io-stdout.ps \
	sparsemat/csrio-stdout.ps \
	sparsemat/csrio-mpiio.ps \
	pneo/pneo_fence.ps \
	pneo/pneo_pscw.ps \
	mutex/mpimutex.ps

#	sparsemat/csrio-pnetcdf.ps
#	life/mlife2d-io-mpiio.ps
#	life/mlife2d-io-pnetcdf.ps

all: $(PSFILES)

# psmerge doesn't work for enscript output, for whatever reason.
# examples.ps: $(PSFILES)
#	psmerge $(PSFILES) > $@

$(PSFILES): %.ps: %.c
	-$(CTOPS) $? $@

clean:
	for dir in $(DIRS) ; do (cd $$dir && make clean) ; done
