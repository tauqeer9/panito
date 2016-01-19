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

#ifndef _GWANI_H_
#define _GWANI_H_

#include "kseq.h"

int is_unknown(char base);
int get_length_of_genome();
int get_number_of_samples();
char ** get_sequence_names();
void check_input_file_and_calc_dimensions(char filename[]);
void calculate_and_output_gwani(char filename[]);
void print_header();
void calc_gwani_between_a_sample_and_everything_afterwards(char filename[], int comparison_index, double * similarity_percentage);

#define DEFAULT_NUM_SAMPLES 65536
#endif
