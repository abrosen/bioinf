# Import smtplib for the actual sending function
import smtplib
# Import the email modules we'll need
from email.mime.text import MIMEText


# station is our server we are connecting to
# it is an SMTP object
# our camels will leave from here and travel with a package
# marked for recipients
def send_camels(station, sender, recipients, package):
    station.sendmail(sender, recipients, package)

MAIL_SERVER = "grid.cs.gsu.edu"


def doTheThing():
    target  = smtplib.SMTP('localhost',1025)
    sender = 'spiderman@friendlyneighborhood.com'
    recv = 'rosen@cs.gsu.edu'
    """
    package = MIMEText("foo this is a test message") 
    package['Subject'] = 'This is a test'
    package['From'] = sender
    package['To'] = recv
    """

    send_camels(target, sender, [recv], '\nhelo')
    target.quit()



if __name__ == '__main__':
    doTheThing()