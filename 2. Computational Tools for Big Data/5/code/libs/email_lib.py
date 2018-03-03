import smtplib
import os
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.MIMEImage import MIMEImage

def send_email(user, pwd, recipient, subject, body):

    gmail_user = user
    gmail_pwd = pwd
    
    FROM = gmail_user
    TO = recipient
    
    # Create the root message and fill in the from, to, and subject headers
    msgRoot = MIMEMultipart('related')
    msgRoot['Subject'] = subject
    msgRoot['From'] = FROM
    msgRoot['To'] = TO
    msgRoot.preamble = 'This is a multi-part message in MIME format.'
    
    # Encapsulate the plain and HTML versions of the message body in an
    # 'alternative' part, so message agents can decide which they want to display.
    msgAlternative = MIMEMultipart('alternative')
    msgRoot.attach(msgAlternative)
    
    msgText = MIMEText('This is the alternative plain text message.')
    msgAlternative.attach(msgText)
    
    # We reference the image in the IMG SRC attribute by the ID we give it below
    msgText = MIMEText('<b>Some <i>HTML</i> text</b> and an image.<br><img src="cid:image1"><br>Nifty!', 'html')
    msgAlternative.attach(msgText)
    
    # This example assumes the image is in the current directory
    fp = open('./asd.png', 'rb')
    msgImage = MIMEImage(fp.read())
    fp.close()
    
    # Define the image's ID as referenced above
    msgImage.add_header('Content-ID', '<image1>')
    msgRoot.attach(msgImage)


    # Prepare actual message

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.ehlo()
        server.starttls()
        server.login(gmail_user, gmail_pwd)
        server.sendmail(FROM, TO, msgRoot.as_string())
        server.close()
        print 'successfully sent the mail'
    except smtplib.SMTPAuthenticationError:
        print "failed to send mail"
        print smtplib.SMTPAuthenticationError

def get_image_as_str(img_dir):
    img_data = open(img_dir, 'rb').read()
    image = MIMEImage(img_data, name =os.path.basename(img_dir))
    return image.as_string()
