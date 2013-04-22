all: FP
	@echo ""
	@echo "Build completed."
	@echo "If you have Sphinx installed you may also want to"
	@echo "build the html documentation by typing 'make documentation'."

FP:
	@echo "Building Fabry-Perot..."
	cd saltfp/fortranfp/; make


bvit:  extension_modules ui_modules
	@echo "Building BVIT modules"

extension_modules:
	@echo "Building extension modules..."
	cd bvittools/accelerate; make

ui_modules:
	@echo "Building graphical user interface modules..."
	cd lib/ui; make

documentation:
	@echo "Building html documentation..."
	cd doc; make html

clean:
	@echo "Removing build ui modules..."
	-cd lib/ui; make clean
	@echo "Removing build extension modules..."
	-cd lib/accelerate; make clean
	@echo "Removing build documentation..."
	-cd doc; make clean
	@echo "Removing byte compiled Python modules..."
	-find . -iname "*.pyc" -exec rm {} \;

