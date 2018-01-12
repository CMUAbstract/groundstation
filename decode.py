#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Top Block
# Generated: Thu May 18 18:03:11 2017
##################################################

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import sprite


class top_block(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Top Block")

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 256000
        self.chip_rate = chip_rate = 48e3

        ##################################################
        # Blocks
        ##################################################
        self.sprite_soft_decoder_c_0 = sprite.soft_decoder_c(.85)
        self.sprite_peak_decimator_ff_0_0 = sprite.peak_decimator_ff(512)
        self.sprite_peak_decimator_ff_0 = sprite.peak_decimator_ff(512)
        self.sprite_correlator_cf_0 = sprite.correlator_cf(266, 267)
        self.blocks_float_to_complex_0 = blocks.float_to_complex(1)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_gr_complex*1, '/tmp/sprite-stream-pipe', False)
        self.blocks_delay_0 = blocks.delay(gr.sizeof_float*1, 256)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_delay_0, 0), (self.sprite_peak_decimator_ff_0_0, 0))    
        self.connect((self.blocks_file_source_0, 0), (self.sprite_correlator_cf_0, 0))    
        self.connect((self.blocks_float_to_complex_0, 0), (self.sprite_soft_decoder_c_0, 0))    
        self.connect((self.sprite_correlator_cf_0, 0), (self.blocks_delay_0, 0))    
        self.connect((self.sprite_correlator_cf_0, 0), (self.sprite_peak_decimator_ff_0, 0))    
        self.connect((self.sprite_peak_decimator_ff_0, 0), (self.blocks_float_to_complex_0, 0))    
        self.connect((self.sprite_peak_decimator_ff_0_0, 0), (self.blocks_float_to_complex_0, 1))    

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate

    def get_chip_rate(self):
        return self.chip_rate

    def set_chip_rate(self, chip_rate):
        self.chip_rate = chip_rate


def main(top_block_cls=top_block, options=None):

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
