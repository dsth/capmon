// #include <sys/types.h>
// #include <sys/socket.h>
// #include <netinet/in.h>
#include <netdb.h>
#include <cstring>
#include <cstdio>
#include <unistd.h> 
#include "email.h"

inline void sendCstring2socket(const char *s, int sock) {

    /*
    bored of ignoring return warning...
    ssize_t ignore;
    ignore = write(sock,s,strlen(s));
    ignore = write(1,s,strlen(s)); 
    */

    write(sock,s,strlen(s));
    write(1,s,strlen(s));
}

inline void readCstringFromSocket(int sock, char *buf, int len) {

    len = read(sock,buf,BUFSIZ);
    write(1,buf,len);

}

struct sockaddr_in server;
struct hostent *hp, *gethostbyname();

bool mailer::send_email(const char * subject, const char * body, const char * _rcpt_to) {

    //y tired of it polluting stdout - this has got to be wrong but can't be bothered to look at it...
    freopen("/dev/null", "w", stdout );

    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock==-1) return false;

    server.sin_family = AF_INET;
    hp = gethostbyname(_smtp);
    if (hp==(struct hostent *) 0) {
        fprintf(stderr, "%s: unknown host\n", _smtp);
        return false;
    }

    memcpy((char *) &server.sin_addr, (char *) hp->h_addr, hp->h_length);

    server.sin_port=htons(25);
    
    if (connect(sock, (struct sockaddr *) &server, sizeof server)==-1) return false;

    readCstringFromSocket(sock,buf,len); 
    sendCstring2socket(HELO,sock); 
    
    readCstringFromSocket(sock,buf,len); 
    sendCstring2socket("MAIL FROM: ",sock); 
    sendCstring2socket(_mail_from,sock);
    sendCstring2socket("\r\n",sock);
    
    readCstringFromSocket(sock,buf,len); 
    sendCstring2socket("VRFY ",sock);
    sendCstring2socket(_mail_from,sock);
    sendCstring2socket("\r\n",sock);
    
    readCstringFromSocket(sock,buf,len); 
    sendCstring2socket("RCPT TO: ",sock); 
    sendCstring2socket(_rcpt_to,sock);
    sendCstring2socket("\r\n",sock);
    
    readCstringFromSocket(sock,buf,len); 
    sendCstring2socket(DATA,sock);
    sendCstring2socket("Subject: ",sock);
    sendCstring2socket(subject,sock);
    sendCstring2socket("\r\n",sock);
    
    readCstringFromSocket(sock,buf,len);
    sendCstring2socket("Mime-Version: 1.0;\r\n",sock);
    sendCstring2socket("Content-Type: text/html; charset=\"ISO-8859-1\";\r\n",sock);
    sendCstring2socket("Content-Transfer-Encoding: 7bit;\r\n",sock);
    sendCstring2socket("<html>\r\n",sock);
    sendCstring2socket("<body>\r\n",sock);
    sendCstring2socket(body,sock);
    sendCstring2socket("</body>\r\n",sock);
    sendCstring2socket("</html>\r\n",sock);
    sendCstring2socket("\r\n",sock);
    sendCstring2socket(".\r\n",sock);
    
    readCstringFromSocket(sock,buf,len); 
    sendCstring2socket(QUIT,sock); 
    
    readCstringFromSocket(sock,buf,len); 

    close(sock);

    freopen("/dev/tty", "w", stdout ); //y just for good measure...

    return true;
}

