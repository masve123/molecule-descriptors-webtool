# The config.py file is where you store configuration settings for your Flask app. This can include:

#Development Configurations: Such as debug mode, local database URIs, etc.
#Production Configurations: Settings for deployment, like production database URIs, secret keys for sessions, etc.
#Testing Configurations: If you plan to write tests for your app.

# En basic versjon, kan måtte endres
import os

class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'
    DEBUG = False

class DevelopmentConfig(Config):
    DEBUG = True

class ProductionConfig(Config):
    DEBUG = False

class TestingConfig(Config):
    TESTING = True
    DEBUG = True
