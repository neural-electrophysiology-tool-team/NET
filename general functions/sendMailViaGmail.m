function sendMailViaGmail(recipientMail,subject,message)

% Modify these two lines to reflect your account and password.
NSKToolBoxMainDir=fileparts(which('identifierOfMainDir4NSKToolBox'));
if exist([NSKToolBoxMainDir filesep 'PCspecificFiles' filesep 'emailsCredentials.txt'],'file')
    emailCredentialsFile = fileread([NSKToolBoxMainDir filesep 'PCspecificFiles' filesep 'emailsCredentials.txt']);
    emailCredentials = regexp(emailCredentialsFile, '\r\n|\r|\n', 'split');
    myaddress = emailCredentials{1};
    mypassword = emailCredentials{2};
else
    message.server.mark@gmail.com
end


setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%{
NET.addAssembly('System.Net')
import System.Net.Mail.*;

mySmtpClient = SmtpClient('smtp.gmail.com');

mySmtpClient.UseDefaultCredentials = false;
mySmtpClient.Credentials = System.Net.NetworkCredential(myaddress, mypassword);

from = MailAddress('from.name@address.com');
to = MailAddress('to.form@address.com');

myMail = MailMessage(from, to);

myMail.Subject = ['Test message: ' datestr(now)];
myMail.SubjectEncoding = System.Text.Encoding.UTF8;
myMail.Body = '<b>Test Mail</b><br>using <b>HTML</b>';
myMail.BodyEncoding = System.Text.Encoding.UTF8;
myMail.IsBodyHtml = true;
mySmtpClient.Send(myMail);
%}

if nargin==3
    try
        sendmail(recipientMail,subject,message);
    catch errorMsg
        disp(['Email was not successfully sent: ' errorMsg.message]);
    end
elseif nargin==2
    try
        sendmail(recipientMail,subject)
    catch errorMsg
        disp(['Email was not successfully sent: ' errorMsg.message]);
    end
else
    error('Either 2 or 3 arguments are required for running sendMailViaGmail function');
end