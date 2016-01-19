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
#include <check.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include "check-panito.h"
#include "gwani.h"

START_TEST (check_small_valid_file_all_same)
{
  check_input_file_and_calc_dimensions("../tests/data/small_all_same.aln");
  fail_unless(get_length_of_genome() == 4);
  fail_unless(get_number_of_samples() == 5);
}
END_TEST

Suite * panito_suite (void)
{
  Suite *s = suite_create ("Creating_panito");

  TCase *tc_gwani = tcase_create ("gwani");
  tcase_add_test (tc_gwani, check_small_valid_file_all_same);
  suite_add_tcase (s, tc_gwani);
  return s;
}

