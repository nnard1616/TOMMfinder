#!/usr/local/bin/python
import os
import pexpect

child = pexpect.spawn('ncftp ftp.ncbi.nlm.nih.gov')
child.expect ('ncftp / > ')
child.sendline ('ls')
