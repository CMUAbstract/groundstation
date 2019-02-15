/* -*- c++ -*- */
/* 
 * Copyright 2014 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "correlator_cf_impl.h"
#include <gnuradio/gr_complex.h>
#include <gnuradio/fft/fft.h>
#include <complex>
#include <cmath>
#include <omp.h>

using namespace std;

namespace gr {
  namespace sprite {

    correlator_cf::sptr
    correlator_cf::make(int prn_id0, int prn_id1)
    {
      return gnuradio::get_initial_sptr
        (new correlator_cf_impl(prn_id0, prn_id1));
    }

    /*
     * The private constructor
     */
    correlator_cf_impl::correlator_cf_impl(int prn_id0, int prn_id1)
      : gr::sync_block("correlator_cf",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(float)))
    {
      set_history(SPRITE_PRN_LENGTH);

      generate_prns(prn_id0, prn_id1);
  
      cc430_modulator(m_prn0, m_template0);
      cc430_modulator(m_prn1, m_template1);
      for (int k = 0; k < SPRITE_PRN_LENGTH; k++)
      {
        m_template0[k] = conj(m_template0[k]);
        m_template1[k] = conj(m_template1[k]);
      }
      
#pragma omp parallel for
      for (int thid = 0; thid < NUM_THREADS; ++thid) {
        m_fft0[thid] = new fft::fft_complex(SPRITE_PRN_LENGTH, true, 1);
        m_fft_buffer_in0[thid] = m_fft0[thid]->get_inbuf();
        m_fft_buffer_out0[thid] = m_fft0[thid]->get_outbuf();

        m_fft1[thid] = new fft::fft_complex(SPRITE_PRN_LENGTH, true, 1);
        m_fft_buffer_in1[thid] = m_fft1[thid]->get_inbuf();
        m_fft_buffer_out1[thid] = m_fft1[thid]->get_outbuf();
      }
    }

    /*
     * Our virtual destructor.
     */
    correlator_cf_impl::~correlator_cf_impl()
    {
    }


    void correlator_cf_impl::cc430_modulator(int* prnBits, gr_complex* baseBand)
    {
      float* diffs = m_buffer_real1;
      float* iBB = m_buffer_real2;
      float* qBB = m_buffer_real3;
      
      //Differentially encode with +/-1 values
      diffs[0] = -2*prnBits[0] + 1;
      for (int k = 1; k < SPRITE_PRN_LENGTH; k++)
      {
        char diff = prnBits[k]-prnBits[k-1];
        if(diff == 0)
        {
          diffs[k] = 1;
        }
        else
        {
          diffs[k] = -1;
        }
      }
      
      //Initialize with offset between I and Q
      iBB[0] = 1;
      qBB[0] = diffs[0];
      qBB[1] = diffs[0];
      
      for(int k = 1; k < SPRITE_PRN_LENGTH-2; k+=2)
      {
        iBB[k] = diffs[k]*iBB[k-1];
        iBB[k+1] = iBB[k];
      }
      iBB[SPRITE_PRN_LENGTH-1] = diffs[SPRITE_PRN_LENGTH-1]*iBB[SPRITE_PRN_LENGTH-2];
      
      for(int k = 2; k < SPRITE_PRN_LENGTH; k+=2)
      {
        qBB[k] = diffs[k]*qBB[k-1];
        qBB[k+1] = qBB[k];
      }
      
      for(int k = 0; k < SPRITE_PRN_LENGTH; k++)
      {
        baseBand[k] = gr_complex(iBB[k]*cos(M_PI/2*k), qBB[k]*sin(M_PI/2*k));
      }
    }

    void correlator_cf_impl::generate_prns(int prn_id0, int prn_id1)
    {
      if(prn_id0 == -2)
      { 
        //Deep copy M-sequence
        for (int k = 0; k < M_SEQUENCE_LENGTH; k++)
        {
          m_prn0[k] = mseq1[k];
        }
      }
      else if(prn_id0 == -1)
      { 
        //Deep copy M-sequence
        for (int k = 0; k < M_SEQUENCE_LENGTH; k++)
        {
          m_prn0[k] = mseq2[k];
        }
      }
      else //if(prn_id >= 0 && prn_id < M_SEQUENCE_LENGTH)
      { 
        //Generate Gold Codes by xor'ing 2 M-sequences in different phases
        for (int k = 0; k < M_SEQUENCE_LENGTH-prn_id0; k++)
        {
          m_prn0[k] = mseq1[k] ^ mseq2[k+prn_id0];
        }
        for (int k = M_SEQUENCE_LENGTH-prn_id0; k < M_SEQUENCE_LENGTH; k++)
        {
          m_prn0[k] = mseq1[k] ^ mseq2[k-M_SEQUENCE_LENGTH+prn_id0];
        }
      }

      m_prn0[SPRITE_PRN_LENGTH-1] = 0; //To pad out the last byte, add a zero to the end

      if(prn_id1 == -2)
      { 
        //Deep copy M-sequence
        for (int k = 0; k < M_SEQUENCE_LENGTH; k++)
        {
          m_prn1[k] = mseq1[k];
        }
      }
      else if(prn_id1 == -1)
      { 
        //Deep copy M-sequence
        for (int k = 0; k < M_SEQUENCE_LENGTH; k++)
        {
          m_prn1[k] = mseq2[k];
        }
      }
      else //if(prn_id >= 0 && prn_id < M_SEQUENCE_LENGTH)
      { 
        //Generate Gold Codes by xor'ing 2 M-sequences in different phases
        for (int k = 0; k < M_SEQUENCE_LENGTH-prn_id1; k++)
        {
          m_prn1[k] = mseq1[k] ^ mseq2[k+prn_id1];
        }
        for (int k = M_SEQUENCE_LENGTH-prn_id1; k < M_SEQUENCE_LENGTH; k++)
        {
          m_prn1[k] = mseq1[k] ^ mseq2[k-M_SEQUENCE_LENGTH+prn_id1];
        }
      }

      m_prn1[SPRITE_PRN_LENGTH-1] = 0; //To pad out the last byte, add a zero to the end

    }

    int
    correlator_cf_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        float *out = (float *) output_items[0];

        // Do <+signal processing+>
#pragma omp parallel for num_threads(NUM_THREADS)
        for(int k = 0; k < noutput_items; ++k) {
          
#if 0
          //Pointwise multiply by baseband template and copy to fft input
          for (int j = 0; j < SPRITE_PRN_LENGTH; ++j)
          {
            m_fft_buffer_in0[j] = m_template0[j]*in[j+k];
            m_fft_buffer_in1[j] = m_template1[j]*in[j+k];
          }
          
          //Take FFT
          m_fft0->execute();
          m_fft1->execute();
#else
          int thid = omp_get_thread_num();
          fft::fft_complex* fft_plans[2] = { m_fft0[thid], m_fft1[thid] };
          gr_complex* in_bufs[2] = { m_fft_buffer_in0[thid], m_fft_buffer_in1[thid] };
          gr_complex* out_bufs[2] = { m_fft_buffer_out0[thid], m_fft_buffer_out1[thid] };
          gr_complex* templates[2] = { m_template0, m_template1 };
          float both_maxmag[2];

#pragma omp parallel for num_threads(2)
          for (int p = 0; p < 2; ++p) {
            //Pointwise multiply by baseband template and copy to fft input
            for (int j = 0; j < SPRITE_PRN_LENGTH; ++j)
            {
              in_bufs[p][j] = templates[p][j]*in[j+k];
            }
            fft_plans[p]->execute();

            float mag = real(out_bufs[p][0]*conj(out_bufs[p][0]));
            float maxmag = mag;
            for (int j = 1; j < SPRITE_PRN_LENGTH; ++j)
            {
              mag = real(out_bufs[p][j]*conj(out_bufs[p][j]));
              if (mag > maxmag)
              {
                maxmag = mag;
              }
            }
            both_maxmag[p] = maxmag;
          }
#endif
          
#if 0
          //Find largest value in FFT
          float mag0 = real(m_fft_buffer_out0[0]*conj(m_fft_buffer_out0[0]));
          float max0 = mag0;
          float index0 = 0;
          for (int j = 1; j < SPRITE_PRN_LENGTH; ++j)
          {
            mag0 = real(out_bufs[0][j]*conj(out_bufs[0][j]));
            if (mag0 > max0)
            {
              max0 = mag0;
              index0 = j;
            }
          }
          float mag1 = real(out_bufs[1][0]*conj(out_bufs[1][0]));
          float max1 = mag1;
          float index1 = 0;
          for (int j = 1; j < SPRITE_PRN_LENGTH; ++j)
          {
            mag1 = real(out_bufs[1][j]*conj(out_bufs[1][j]));
            if (mag1 > max1)
            {
              max1 = mag1;
              index1 = j;
            }
          }
          out[k] = max1 >= max0 ? sqrt(max1) : -sqrt(max0);
#else
          out[k] = both_maxmag[1] >= both_maxmag[0] ? sqrt(both_maxmag[1]) : -sqrt(both_maxmag[0]);
#endif
          
        }

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace sprite */
} /* namespace gr */

