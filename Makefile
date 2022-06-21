all: install

format:
	clang-format -i *.hpp

install: 
	rm -rf ~/toofus
	mkdir ~/toofus
	cp ./*.hpp ~/toofus
	
	