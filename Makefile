DESTDIR ?= /

SCRIPTS = \
	SpriteReceiver2.py \
	SpriteRecorder_RTL.py \
	SpriteDecoder.py \
	WavToRaw.py \

all : $(SCRIPTS)

# NOTE: this assumes the .grc file's 'id' property is set to the filename
%.py : %.grc
	grcc -d . $^

# This is hacky.. perhaps a better way would be to create a 'spritectl' python package
# and install it using setuptools
 
install : $(SCRIPTS)
	for s in $(SCRIPTS); do install -Dm755 $$s $(DESTDIR)/usr/bin/$$s; done

clean :
	rm -f $(SCRIPTS)
