#!/usr/bin/python

# Import the CGI module
import cgi

# Required header that tells the browser how to render the HTML.
print "Content-Type: text/html\n\n"

# Define function to generate HTML form.


def generate_form():
    print "<HTML>\n"
    print "<HEAD>\n"
    print "\t<TITLE>Info Form</TITLE>\n"
    print "</HEAD>\n"
    print "<BODY BGCOLOR = white>\n"
    print "\t<H3>Please, enter your name and age.</H3>\n"
    print "\t<TABLE BORDER = 0>\n"
    print "\t\t<FORM METHOD = post ACTION = \
    \"example_7.1.cgi\">\n"
    print "\t\t<TR><TH>Name:</TH><TD><INPUT type = text \
    name = \"name\"></TD><TR>\n"
    print "\t\t<TR><TH>Age:</TH><TD><INPUT type = text name = \
    \"age\"></TD></TR>\n"
    print "\t</TABLE>\n"
    print "\t<INPUT TYPE = hidden NAME = \"action\" VALUE = \
    \"display\">\n"
    print "\t<INPUT TYPE = submit VALUE = \"Enter\">\n"
    print "\t</FORM>\n"
    print "</BODY>\n"
    print "</HTML>\n"

# Define function display data.


def display_data(name, age):
    print "<HTML>\n"
    print "<HEAD>\n"
    print "\t<TITLE>Info Form</TITLE>\n"
    print "</HEAD>\n"
    print "<BODY BGCOLOR = white>\n"
    print name, ", you are", age, "years old."
    print "</BODY>\n"
    print "</HTML>\n"

# Define main function.


def main():
    form = cgi.FieldStorage()
    if ("action" in form and "name" in form and "age" in form):
        if (form["action"].value == "display"):
            display_data(form["name"].value, form["age"].value)
        else:
            generate_form()

# Call main function.
main()
