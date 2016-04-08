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

void fast_calculate_gwani(char filename[])
{
  check_input_file_and_calc_dimensions(filename);
  print_header();
  
  // initialise space to store entire genome
  int i;
  int j;
  char ** comparison_sequence;
  comparison_sequence = calloc(get_number_of_samples() + 1, sizeof(char *));
  for(i=0; i < get_number_of_samples(); i++)
  {
    comparison_sequence[i] = calloc(get_length_of_genome() + 1, sizeof(char));
  }
  
  // Store all sequences in a giant array - eek
  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(filename, "r");
  seq = kseq_init(fp);
  i =0;
  while ((l = kseq_read(seq)) >= 0) {
    memcpy(comparison_sequence[i], seq->seq.s, strlen(seq->seq.s)+1);
    i++;
  }
  
  for(j = 0; j< number_of_samples; j++)
  {
    for(i = 0; i < length_of_genome; i++)
    {
      //standardise the input so that case doesnt matter and replace unknowns with single type
      comparison_sequence[j][i] = toupper(comparison_sequence[j][i]);
      if(is_unknown(comparison_sequence[j][i]))
      {
        comparison_sequence[j][i] = 'N';
      }
    }
  }

  for(i = 0; i < number_of_samples; i++)
  {
    printf("%s",sequence_names[i]);
    double * similarity_percentage;
    similarity_percentage = calloc(number_of_samples + 1 , sizeof(double));

    calc_gwani_between_a_sample_and_everything_afterwards_memory(comparison_sequence, i ,similarity_percentage);
    
    for(j = 0; j < number_of_samples; j++)
    {
      if(similarity_percentage[j] < 0)
      {
        printf("\t-");
      }
      else
      {
        printf("\t%f",similarity_percentage[j]);
      }
    }
    printf("\n");
    free(similarity_percentage);
  }
  
}
  


void calc_gwani_between_a_sample_and_everything_afterwards_memory(char ** comparison_sequence,int comparison_index, double * similarity_percentage)
{
  int current_index = 0;
  int i;
  int l;
  int j;
  int bases_in_common;
  int length_without_gaps;

  for(j = 0; j< number_of_samples; j++)
  {
    if(current_index < comparison_index)
    {
      similarity_percentage[current_index] = -1;
    }
    else if(current_index == comparison_index)
    {
      similarity_percentage[current_index] = 100;
    }
    else
    {
      bases_in_common = 0;
      length_without_gaps =length_of_genome;
      for(i = 0; i < length_of_genome; i++)
      {
        
        if(comparison_sequence[comparison_index][i] == 'N' || comparison_sequence[j][i] == 'N' )
        {
          length_without_gaps--;
        } 
        else if(comparison_sequence[comparison_index][i] == comparison_sequence[j][i] )
        {
          bases_in_common++;
        }
      }
      if(length_without_gaps > 0)
      {
          similarity_percentage[current_index] = (bases_in_common*100.0)/length_without_gaps;
      } 
      else
      {
        similarity_percentage[current_index] = 0;
      }
      
    }
    current_index++;
  }
}


void calculate_and_output_gwani(char filename[])
{
  check_input_file_and_calc_dimensions(filename);
  print_header();
  
  int i;
  int j;
  for(i = 0; i < number_of_samples; i++)
  {
    printf("%s",sequence_names[i]);
    double * similarity_percentage;
    similarity_percentage = calloc(number_of_samples + 1 , sizeof(double));
    calc_gwani_between_a_sample_and_everything_afterwards(filename, i, similarity_percentage);
    
    for(j = 0; j < number_of_samples; j++)
    {
      if(similarity_percentage[j] < 0)
      {
        printf("\t-");
      }
      else
      {
        printf("\t%f",similarity_percentage[j]);
      }
    }
    printf("\n");
    free(similarity_percentage);
  }
}

void print_header()
{
  int i;
  for(i = 0; i < number_of_samples; i++)
  {
    printf("\t%s",sequence_names[i]);
  }
  printf("\n");
}

void calc_gwani_between_a_sample_and_everything_afterwards(char filename[],int comparison_index, double * similarity_percentage)
{
  int current_index = 0;
  int i;
  int l;
  int bases_in_common;
  int length_without_gaps;
  gzFile fp;
  kseq_t *seq;

  char * comparison_sequence;
  comparison_sequence = calloc(length_of_genome + 1, sizeof(char));
  
  fp = gzopen(filename, "r");
  seq = kseq_init(fp);
  
  while ((l = kseq_read(seq)) >= 0) {
    if(current_index < comparison_index)
    {
      similarity_percentage[current_index] = -1;
    }
    else if(current_index == comparison_index)
    {
      similarity_percentage[current_index] = 100;
      for(i = 0; i < length_of_genome; i++)
      {
          //standardise the input so that case doesnt matter and replace unknowns with single type
          comparison_sequence[i] = toupper(seq->seq.s[i]);
          if(is_unknown(comparison_sequence[i]))
          {
            comparison_sequence[i] = 'N';
          }
      }
    }
    else
    {
      bases_in_common = 0;
      length_without_gaps =length_of_genome;
      for(i = 0; i < length_of_genome; i++)
      {
        
        if(comparison_sequence[i] == 'N' || is_unknown(seq->seq.s[i]) )
        {
          length_without_gaps--;
        } 
        else if(comparison_sequence[i] == toupper(seq->seq.s[i]) && ! is_unknown(seq->seq.s[i]))
        {
          bases_in_common++;
        }
      }
      if(length_without_gaps > 0)
      {
          similarity_percentage[current_index] = (bases_in_common*100.0)/length_without_gaps;
      } 
      else
      {
        similarity_percentage[current_index] = 0;
      }
      
    }
    current_index++;
  }
  
  kseq_destroy(seq);
  gzclose(fp);
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

//void cleanup_memory
//{
//  free(sequence_names);
//}

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
