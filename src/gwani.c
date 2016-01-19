/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2016  Wellcome Trust Sanger Institute
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 3
 *  of the License, or (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <regex.h>
#include <sys/types.h>
#include "kseq.h"
#include "gwani.h"

KSEQ_INIT(gzFile, gzread)

int length_of_genome;
int number_of_samples;
char ** sequence_names;


int get_length_of_genome()
{
    return length_of_genome;
}

int get_number_of_samples()
{
    return number_of_samples;
}

char ** get_sequence_names()
{
    return sequence_names;
}

void check_input_file_and_calc_dimensions(char filename[])
{
  int i;
  int l;
  number_of_samples = 0; 
  length_of_genome = 0;
  gzFile fp;
  kseq_t *seq;
  
  fp = gzopen(filename, "r");
  seq = kseq_init(fp);
  sequence_names = calloc(DEFAULT_NUM_SAMPLES, sizeof(char*));

  // First pass of the file get the length of the alignment, number of samples and sample names
  while ((l = kseq_read(seq)) >= 0) {
    if(number_of_samples == 0)
    {
        length_of_genome = seq->seq.l;
    }
    else if(length_of_genome != seq->seq.l)
    {
      fprintf(stderr, "Alignment %s contains sequences of unequal length. Expected length is %i but got %i in sequence %s\n\n",filename, (int) length_of_genome, (int) seq->seq.l,seq->name.s);
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    // The sample name is initially set to a large number but make sure this can be increased dynamically
   if(number_of_samples >= DEFAULT_NUM_SAMPLES)
   {
     sequence_names = realloc(sequence_names, (number_of_samples + 1) * sizeof(char*));
   }
   sequence_names[number_of_samples] = calloc(FILENAME_MAX,sizeof(char));
   strcpy(sequence_names[number_of_samples], seq->name.s);
   
   number_of_samples++;
  }

  kseq_destroy(seq);
  gzclose(fp);
  return;
}

int is_unknown(char base)
{
  switch (base) {
    case 'N':
    case 'n':
    case '-':
    case '?':
      return 1;
    default:
      return 0;
  }
}
