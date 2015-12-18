from django.conf.urls import patterns, url
from main.views import *

urlpatterns = patterns('',
                       url(r'^$', home, name='home'),
)
