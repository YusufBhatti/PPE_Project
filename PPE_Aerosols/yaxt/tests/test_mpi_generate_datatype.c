/**
 * @file test_mpi_generate_datatype.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"

#undef FIXED_MPI_ON_BLIZZARD

#ifdef FIXED_MPI_ON_BLIZZARD

static int
test_datatype_int(MPI_Datatype datatype, int recv_count, int * send_data,
                  int * ref_recv_data) {

  int recv_data[recv_count];
  MPI_Status    status;


  //datatype
  xt_mpi_call(MPI_Sendrecv(send_data, 1, datatype, 0, 0, recv_data, recv_count,
                           MPI_INT, 0, 0, MPI_COMM_WORLD, &status),
              MPI_COMM_WORLD);

  for (int i = 0; i < recv_count; ++i)
    if (recv_data[i] != ref_recv_data[i])
      return 1;

  return 0;
}

#else

static int
test_datatype_int(MPI_Datatype datatype, int recv_count, int * send_data,
                  int * ref_recv_data) {

  int recv_data[recv_count];
  MPI_Status    status;
  MPI_Request request;

  xt_mpi_call(MPI_Irecv(recv_data, recv_count, MPI_INT, 0, 0, MPI_COMM_WORLD,
                        &request), MPI_COMM_WORLD);

  xt_mpi_call(MPI_Send(send_data, 1, datatype, 0, 0, MPI_COMM_WORLD),
              MPI_COMM_WORLD);

  xt_mpi_call(MPI_Wait(&request, &status), MPI_COMM_WORLD);

  for (int i = 0; i < recv_count; ++i)
    if (recv_data[i] != ref_recv_data[i])
      return 1;

  return 0;
}

#endif

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // tests with empty data

    MPI_Datatype datatype;

    int displacements[1] = {0};
    int blocklengths[1] = {0};
    int count = 0;
    MPI_Datatype old_type = MPI_INT;

    datatype = xt_mpi_generate_datatype(displacements, count, old_type,
                                        MPI_COMM_WORLD);

    if (datatype != MPI_DATATYPE_NULL)
      PUT_ERR("error in xt_mpi_generate_datatype (count == 0)\n");

    datatype = xt_mpi_generate_datatype_block(displacements, blocklengths,
                                              count, old_type, MPI_COMM_WORLD);

    if (datatype != MPI_DATATYPE_NULL)
      PUT_ERR("error in xt_mpi_generatxt_mpi_generate_datatype_blocke_datatype"
              " (count == 0)\n");
  }

  { // tests with block_count == 1

    { // blocklength == 1 and displacement = 0
      MPI_Datatype datatype;

      int displacements[1] = {0};
      int blocklengths[1] = {1};
      int count = 1;
      MPI_Datatype old_type = MPI_INT;

      int send_data[3] = {123,234,345};
      int ref_recv_data[1] = {123};

      datatype = xt_mpi_generate_datatype(displacements, count, old_type,
                                          MPI_COMM_WORLD);

      if (test_datatype_int(datatype, count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(displacements, blocklengths,
                                                count, old_type,
                                                MPI_COMM_WORLD);

      if (test_datatype_int(datatype, count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by"
                " xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // blocklength == 1 and displacement = 1
      MPI_Datatype datatype;

      int displacements[1] = {1};
      int blocklengths[1] = {1};
      int count = 1;
      MPI_Datatype old_type = MPI_INT;

      int send_data[3] = {123,234,345};
      int ref_recv_data[1] = {234};

      datatype = xt_mpi_generate_datatype(displacements, count, old_type,
                                          MPI_COMM_WORLD);

      if (test_datatype_int(datatype, count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(displacements, blocklengths,
                                                count, old_type,
                                                MPI_COMM_WORLD);

      if (test_datatype_int(datatype, count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);

    }

    { // blocklength == 2 and displacement = 0
      MPI_Datatype datatype;

      int element_displacements[2] = {0,1};
      int block_displacements[1] = {0};
      int blocklengths[1] = {2};
      int element_count = 2;
      int block_count = 1;
      MPI_Datatype old_type = MPI_INT;

      int send_data[3] = {123,234,345};
      int ref_recv_data[2] = {123, 234};

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // blocklength == 2 and displacement = 1
      MPI_Datatype datatype;

      int element_displacements[2] = {1,2};
      int block_displacements[1] = {1};
      int blocklengths[1] = {2};
      int element_count = 2;
      int block_count = 1;
      MPI_Datatype old_type = MPI_INT;

      int send_data[3] = {123,234,345};
      int ref_recv_data[2] = {234,345};

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths,
                                                block_count, old_type,
                                                MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }
  }

  { // tests with block_count == 3
    { // block length == {1,1,1} and block displacement = {0,2,4}
      MPI_Datatype datatype;

      int element_displacements[3] = {0,2,4};
      int block_displacements[3] = {0,2,4};
      int blocklengths[3] = {1,1,1};
      int element_count = 3;
      int block_count = 3;
      MPI_Datatype old_type = MPI_INT;

      int send_data[5] = {123,234,345,456,567};
      int ref_recv_data[3] = {123,345,567};

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,2,2} and block displacement = {0,3,6}
      MPI_Datatype datatype;

      int element_displacements[6] = {0,1,3,4,6,7};
      int block_displacements[3] = {0,3,6};
      int blocklengths[3] = {2,2,2};
      int element_count = 6;
      int block_count = 3;
      MPI_Datatype old_type = MPI_INT;

      int send_data[8] = {0,1,2,3,4,5,6,7};
      int ref_recv_data[6] = {0,1,3,4,6,7};

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,2,2} and block displacement = {1,4,7}
      MPI_Datatype datatype;

      int element_displacements[6] = {1,2,4,5,7,8};
      int block_displacements[3] = {1,4,7};
      int blocklengths[3] = {2,2,2};
      int element_count = 6;
      int block_count = 3;
      MPI_Datatype old_type = MPI_INT;

      int send_data[9] = {0,1,2,3,4,5,6,7,8};
      int ref_recv_data[6] = {1,2,4,5,7,8};

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,2,2} and block displacement = {7,4,1}
      MPI_Datatype datatype;

      int element_displacements[6] = {7,8,4,5,1,2};
      int block_displacements[3] = {7,4,1};
      int blocklengths[3] = {2,2,2};
      int element_count = 6;
      int block_count = 3;
      MPI_Datatype old_type = MPI_INT;

      int send_data[9] = {0,1,2,3,4,5,6,7,8};
      int ref_recv_data[6] = {7,8,4,5,1,2};

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,2,2} and block displacement = {0,3,7}
      MPI_Datatype datatype;

      int element_displacements[6] = {0,1,3,4,7,8};
      int block_displacements[3] = {0,3,7};
      int blocklengths[3] = {2,2,2};
      int element_count = 6;
      int block_count = 3;
      MPI_Datatype old_type = MPI_INT;

      int send_data[9] = {0,1,2,3,4,5,6,7,8};
      int ref_recv_data[6] = {0,1,3,4,7,8};

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,3,4} and block displacement = {0,3,7}
      MPI_Datatype datatype;

      int element_displacements[9] = {0,1,3,4,5,7,8,9,10};
      int block_displacements[3] = {0,3,7};
      int blocklengths[3] = {2,3,4};
      int element_count = 9;
      int block_count = 3;
      MPI_Datatype old_type = MPI_INT;

      int send_data[11] = {0,1,2,3,4,5,6,7,8,9,10};
      int ref_recv_data[9] = {0,1,3,4,5,7,8,9,10};

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }
  }


  { // check two encoded trivial vectors
    MPI_Datatype datatype;

    int element_displacements[4] = {2,4,7,9};
    int element_count = 4;
    MPI_Datatype old_type = MPI_INT;

    int send_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
    int ref_recv_data[4] = {2,4,7,9};
    datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

    if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
      PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

    MPI_Type_free(&datatype);

  }

  { // check two encoded non-trivial vectors
    MPI_Datatype datatype;

    int element_displacements[8] = {1,4,7,10, 11,14,17,20};
    int element_count = 8;
    MPI_Datatype old_type = MPI_INT;

    int send_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
    int ref_recv_data[8] = {1,4,7,10, 11,14,17,20};
    datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);
    if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
      PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

    MPI_Type_free(&datatype);
  }

  { // check reverse block
    MPI_Datatype datatype;

    int element_displacements[14] = {14,13,12,11,10,9,8,7,6,5,4,3,2,1};

    int element_count = 14;
    MPI_Datatype old_type = MPI_INT;

    int send_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
    int ref_recv_data[14] = {14,13,12,11,10,9,8,7,6,5,4,3,2,1};
    datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);
    if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
      PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

    MPI_Type_free(&datatype);

  }

  xt_finalize();
  xt_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
}