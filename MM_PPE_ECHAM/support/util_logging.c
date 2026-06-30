#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>

#include <sys/socket.h>
#include <arpa/inet.h>

#include <netdb.h>

int sock;                         /* socket descriptor */
struct sockaddr_in ServAddr;      /* server address */

int open_network_logging(char *hostname, int port)
{ 
  unsigned short ServPort;
  struct hostent *host;
  /* TB: 
     FMI Cray XC30 crashes here without __XC30__ */
#if defined (__XT3__) || defined (__XE6__)  || defined (__XC30__)
  host = NULL;
#else
  host = gethostbyname(hostname);
#endif

  if ((host == NULL))
    {
      perror("gethostbyname() failed");
      return(1);
    }

  ServPort = (unsigned int) port;   /* assign port to right data type*/
 
  /* Create a reliable, stream socket using TCP */

  if ((sock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0)
    {
      perror("socket() failed");
      return(1);
    }
  
  /* Construct the server address structure */

  memset(&ServAddr, 0, sizeof(ServAddr));

  ServAddr.sin_family      = AF_INET;
  ServAddr.sin_addr.s_addr = inet_addr(inet_ntoa(*(struct in_addr *) (host->h_addr_list[0])));
  ServAddr.sin_port        = htons(ServPort);
  
  /* Establish the connection to the echo server */

  if (connect(sock, (struct sockaddr *) &ServAddr, sizeof(ServAddr)) < 0)
    {
      perror("connect() failed");
      return(1);
    }

  return(0);
}

int send_network_log_message(char *message)
{
  int MessageLen;
  char *SendBuffer;

  MessageLen = strlen(message);          /* Determine message length */

  SendBuffer = (char *) malloc(sizeof(char)*(MessageLen+3));

  strncpy(SendBuffer, message, (size_t) MessageLen+1);
  strcat(SendBuffer, "\r\n");

  /* Send the string to the server */

  if (send(sock, SendBuffer, MessageLen+2, 0) != MessageLen+2)
    {  
      perror("send() sent a different number of bytes than expected");
      return(1);
    }
  return (0);

}

void close_network_logging(void)
{
  close(sock);
}
