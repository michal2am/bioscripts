import os, time
import smtplib
from pyngrok import exception
from pyngrok.conf import PyngrokConfig
from pyngrok import ngrok, conf
import smtplib, ssl

# ngrok setup

conf.get_default().region = 'eu'

ssh_tunnel = ngrok.connect(22, "tcp", options={'region': 'eu'})
ngrok_process = ngrok.get_ngrok_process()
ngrok_tunnels = ngrok.get_tunnels()

context = ssl.create_default_context()

serw_ngrok = str(ngrok_tunnels[0])[20]
port_ngrok = str(ngrok_tunnels[0])[38:43]

# mailing stuff:

port = 587
smtp_server = "smtp-mail.outlook.com"
password = 'Gaba1111'
from_adr = "tatryworkstation@outlook.com"

message = """\
Subject: NGROK

server:{} port:{}""".format(serw_ngrok, port_ngrok)

send_mail = True

if send_mail:
    with smtplib.SMTP(smtp_server, port) as server:
        server.starttls(context=context)
        server.login(from_adr, password)
        server.sendmail(from_adr, "michaladammichalowski@gmail.com", message)

# ngrok loop

try:
    # Block until CTRL-C or some other terminating event
    ngrok_process.proc.wait()

except KeyboardInterrupt:
    print(" Shutting down server.")
    ngrok.kill()
