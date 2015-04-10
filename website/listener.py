from http.server import HTTPServer, BaseHTTPRequestHandler
import cgi
import time
import json

rna_form_page = b"""

<html>
<header>
<title>siRNA Project</title>
</header>
<body>
<h1>siRNA Project</h1>

<form  enctype="multipart/form-data" method="POST">


<input type="file" name="data"/>
<input type="text" name="email"/>

<input type="submit" value="Submit File"/>
</form>

</body>
</html>

"""

class InitialHandler(BaseHTTPRequestHandler):
    
    def do_GET(self):
        self.send_response(200)
        self.send_header('Content-type','text/html')
        self.end_headers()
        self.wfile.write(rna_form_page)
        
    def do_POST(self):
        form = cgi.FieldStorage(
            fp=self.rfile,
            headers=self.headers,
            environ={'REQUEST_METHOD':'POST','CONTENT_TYPE':self.headers['Content-Type'],})
        filename = form['data'].filename
        data = form['data'].file.read()
        email  = form['email'].file.read()
        print(filename)
        print(data)
        print(email)
        
        
        open("uploads/"+ str(int(time.time()*1000))+filename, "wb").write(data)
        open('emails.txt', "a").write(filename + " " + email + "\n")
        self.send_response(202)
        self.send_header('Content-type','text/html')
        self.end_headers()
        self.wfile.write(b"Success. Go away and stop bothering me and I'll email you later.")


def run(server_class=HTTPServer, handler_class=InitialHandler):
    server_address = ('', 8000)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()


run()
