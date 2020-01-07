import smtplib, ssl

def send_finished_script_email(custom_message = ''):
    port = 465  # For SSL

    with open("email-credentials") as credentials_file:
        sender_email = credentials_file.readline().strip('\n')
        password = credentials_file.readline().strip('\n')
        receiver_email = credentials_file.readline()


    message = """\
From: {}
To: {}
Subject: Script finished

This message is sent from Python.

{}""".format(sender_email, receiver_email,custom_message)

    # Create a secure SSL context
    context = ssl.create_default_context()

    with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
        server.login("m.hynes.servernotifier@gmail.com", password)
        server.sendmail(sender_email,receiver_email, message)
        print("Sent email")

if __name__ == "__main__":
    send_finished_script_email()