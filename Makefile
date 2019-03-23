

exsort: exsort.c exsort.h gene_core.c gene_core.h
	gcc exsort.c gene_core.c -lz -lpthread -lc -o exsort

clean:
	rm -f exsort
