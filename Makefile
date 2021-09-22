all: install

format:
	clang-format -i *.hpp

install: 
	rm -rf /usr/local/include/toofus
	mkdir /usr/local/include/toofus
	cp ./*.hpp /usr/local/include/toofus
	
	