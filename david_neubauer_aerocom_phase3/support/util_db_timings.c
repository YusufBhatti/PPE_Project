#include <stdio.h>
#include <stdlib.h>

#include "config.h"

#ifdef HAVE_POSTGRESQL
#include "libpq-fe.h"
#endif

const int success = 0;
const int failed  = 1;

#if HAVE_POSTGRESQL == 1

static PGconn *connection;

int open_db_timings(char *connectionInfo)
{
  connection = PQconnectdb(connectionInfo);
  if (PQstatus(connection) != CONNECTION_OK)
    {
      fprintf(stderr, "Connection not established: %s\n", 
	      PQerrorMessage(connection));
      PQfinish(connection);
      return failed;
    }
  return success;
}

int begin_transaction_block(void)
{
  PGresult *res;

  res = PQexec(connection, "BEGIN");
  if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
      fprintf(stderr, "BEGIN command failed\n");
      PQclear(res);
      PQfinish(connection); 
      return failed;
    }
  PQclear(res);
  return success;
}

int commit_transaction_block(void)
{
  PGresult *res;

  res = PQexec(connection, "COMMIT");
  if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
      fprintf(stderr, "COMMIT command failed\n");
      PQclear(res);
      PQfinish(connection); 
      return failed;
    }
  PQclear(res);
  return success;
}

int insert_experiment(char *experiment, char *outpath, char *host, char *os, char *user)
{
  PGresult *res;
  char query[1280];

  snprintf(query, 1280, 
	   "INSERT INTO experiments VALUES "
	   "( '%s', '%s', '%s', '%s', '%s' );",
	   experiment, outpath, host, os, user);
#ifdef _DEBUG
  fprintf(stderr, "insert experiment %s\n", query);
#endif
  res = PQexec(connection, query);
  if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
      fprintf(stderr, "Insertion of experiment failed: %s\n", 
	      PQerrorMessage(connection));
      PQclear(res);
      /* PQfinish(connection); */
      return failed;
    }
  PQclear(res);
  return success;
}

int insert_job(char *label, char *experiment, 
	       int nproca, int nprocb, int nprocio, int nthreads, 
	       char *job_run_time_unit, int job_run_time, int timestep,
	       int truncation, int levels, int ncycles)
{
  PGresult *res;
  char query[512];

  snprintf(query, 512, 
	   "INSERT INTO jobs VALUES "
	   "( '%s', '%s', %d, %d, %d, %d, '%s', %d, %d, %d, %d, %d, NOW() );",
	   label, experiment, nproca, nprocb, nprocio, nthreads,
	   job_run_time_unit, job_run_time, timestep, truncation, levels, 
	   ncycles);
#ifdef _DEBUG
  fprintf(stderr, "insert job %s\n", query);
#endif
  res = PQexec(connection, query);
  if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
      fprintf(stderr, "Insertion of job failed: %s\n", 
	      PQerrorMessage(connection));
      PQclear(res);
      /* PQfinish(connection); */
      return failed;
    }
  PQclear(res);
  return success;
}

int insert_timing(char *label, char *job_label,
		  double minval, double avgval, double maxval, double sumval, 
		  double efficiency)
{
  PGresult *res;
  char query[256];

  snprintf(query, 256, 
	   "INSERT INTO timings VALUES "
	   "( '%s', '%s', %lf, %lf, %lf, %lf, %lf );",
	   label, job_label, minval, avgval, maxval, sumval, efficiency);
#ifdef _DEBUG
  fprintf(stderr, "insert timing %s\n", query);  
#endif
  res = PQexec(connection, query);
  if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
      fprintf(stderr, "Insertion of timing record failed: %s\n", 
	      PQerrorMessage(connection));
      PQclear(res);
      /* PQfinish(connection); */
      return failed;
    }
  PQclear(res);
  return success;
}

void close_db_timings(void)
{
  PQfinish(connection);
  return;
}

#else

int open_db_timings(void)
{
  return success;
}

int begin_transaction_block(void)
{
  return success;
}

int commit_transaction_block(void)
{
  return success;
}

int insert_experiment(char *experiment)
{
  return success;
}

int insert_job(char *label, char *experiment, 
	       int nproca, int nprocb, int nprocio, int nthreads, 
	       char *job_run_time_unit, int job_run_time, int timestep,
	       int truncation, int levels, int ncycles)
{
  return success;
}

int insert_timing(char *label, char *job_label,
		  double minval, double avgval, double maxval, double sumval, 
		  double efficiency)
{
  return success;
}

void close_db_timings(void)
{
  return;
}

#endif
