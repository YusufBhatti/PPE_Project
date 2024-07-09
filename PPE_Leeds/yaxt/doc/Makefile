
html/index.html: ../configure.ac Doxyfile.in
	cd .. && ./scripts/generate_doc.pl configure.ac \
		doc/Doxyfile.in doc/Doxyfile && \
		doxygen doc/Doxyfile


clean:
	rm -rf html/

