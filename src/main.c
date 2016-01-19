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
#include <ctype.h>
#include <unistd.h>
#include "config.h"
#include "gwani.h"

#define PROGRAM_NAME "panito"
#define PROGRAM_VERSION PACKAGE_VERSION

static void print_usage()
{
	printf("Usage: panito [-hV] <file>\n");
	printf("This program calculates the genome wide ANI for a multiFASTA alignment.\n");
	printf(" -h		this help message\n");
	printf(" -V		print version and exit\n");
	printf(" <file>		input alignment file which can optionally be gzipped\n");
}

static void print_version()
{
  printf("%s %s\n", PROGRAM_NAME, PROGRAM_VERSION);
}

int main (int argc, char **argv) {
	char multi_fasta_filename[FILENAME_MAX] = {""};
	char output_filename[FILENAME_MAX] = {""};

	int c;
	int index;
  int output_multi_fasta_file = 0;
	
	 while ((c = getopt (argc, argv, "ho:V")) != -1)
      switch (c)
        {
        case 'V':
          print_version();
          return 0;
	      case 'o':
          strncpy(output_filename, optarg, FILENAME_MAX);
	        break;
        case 'h':
          print_usage();
          return 0;
      default:
        output_multi_fasta_file = 1;
      }
  
  if(optind < argc)
  {
    strncpy(multi_fasta_filename, argv[optind], FILENAME_MAX); 
    calculate_and_output_gwani(multi_fasta_filename);
  }
  else
  {
    print_usage();
  }
    
	return 0;
}


