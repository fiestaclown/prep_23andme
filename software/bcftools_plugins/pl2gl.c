/*  pl2gl.c converts FORMAT/PL to FORMAT/GL. removes FORMAT/PL

Copyright (C) 2015 Illumina

Author: Jared O'Connell <joconnell@illumina.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <math.h>

bcf_hdr_t *in_hdr,*out_hdr;
float *gl = NULL;int ngl;
int32_t *pl = NULL;int npl;
int n;
const char *about(void)
{
  return "converts PL to GL\n";
}



int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{

  in_hdr  = in;
  out_hdr = out;
  n =  bcf_hdr_nsamples(in_hdr);
  ngl=3*n;
  npl=ngl;
  gl = (float *)malloc(ngl*sizeof(float));
  pl = (int32_t *)malloc(npl*sizeof(int32_t));
  bcf_hdr_append(out_hdr,"##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"three log10-scaled likelihoods for RR,RA,AA genotypes\">");
  return 0;
}


bcf1_t *process(bcf1_t *rec)
{
  int i;
  if(rec->n_allele==2) {
    assert( bcf_get_format_int32(in_hdr, rec, "PL", &pl, &npl) == (n*3));
    for(i=0;i<ngl;i++) {
      gl[i]=0;
      if(pl[i]>0) {
	gl[i] = (float)pl[i];
	gl[i] /= -10.;
      }
    }
    bcf_update_format_float(out_hdr,rec,"GL",gl,ngl);    
    bcf_update_format_float(out_hdr,rec,"PL",NULL,0);
  }
  return rec;  
}

void destroy(void)
{
  free(gl);
  free(pl);
}
