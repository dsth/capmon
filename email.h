#ifndef EMAIL_H 
#define EMAIL_H

#define HELO "HELO 127.0.0.1\r\n"
#define DATA "DATA\r\n"
#define QUIT "QUIT\r\n"

class mailer {
    const char *_smtp;
    const char *_mail_from;

    int sock; 
    char buf[BUFSIZ+1];
    int len;

  public:
    mailer(const char *smtp, const char *from) : _smtp(smtp), _mail_from(from) {};
    bool send_email(const char *, const char *, const char *);
};

#endif           

