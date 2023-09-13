# molecule-descriptors-webtool
A simple, user-friendly website where users can input molecules (a string) and use RDkit, an open source Python package, to calculate various descriptors.

## Oppdateres etter hvert
- For å hente nyeste versjon kjøres ```git pull``` i terminalen
- I den første versjonen har jeg laget et utkast til backend-delen av prosjektet. ```requirements.txt``` inneholder informasjon om de riktige pakkene som trengs (foreløpig har jeg kun lagt inn RDkit og Flask. For å laste ned de riktige pakkene kan man kjøre ```pip install -r requirements.txt``` i terminalen.

# Introduction
A simple, user-friendly website where users can input molecules (a string) and use RDkit, an open source Python package, to calculate various descriptors.


## For å hente nyeste versjon kjøres git pull i terminalen
I den første versjonen har jeg laget et utkast til backend-delen av prosjektet. For å laste ned de riktige pakkene kan man kjøre pip install -r requirements.txt i terminalen etter å ha pullet prosjektet.



- /static: 
This directory is used to store static files, such as CSS files, JavaScript files, and images. These files don't change and are served as-is to the client.
For example, if you have a CSS file to style your web application, you would place it in the static directory.
/templates:

- /templates
This directory is used to store HTML templates. Flask uses a templating engine called Jinja2 that allows you to embed Python code within HTML files. This is how dynamic content is rendered.
For example, if you want to display the molecular descriptors in an HTML page, you'd have an HTML template in the templates directory that Flask fills with the appropriate data before sending it to the client.

- The config.py file is where you store configuration settings for your Flask app. This can include:

- requirements.txt is the file that holds the dependencies (packages) that are needed. This includes RDkit and Flask.
