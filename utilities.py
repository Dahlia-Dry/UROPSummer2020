import os
import subprocess
from paramiko import SSHClient
from scp import SCPClient
import scp
from astropy.io import fits
import glob
import shutil
import pandas as pd
import numpy as np
import sys

def satisfies_constraints(constraints):
    """ Returns list of #### strings representing the file numbers that satisfy
    the given dict of constraints.
    ...
    Opens an SSH connection to the mitastrometry server, grabs all the Rb and
    .LOG files, shoves them into temp directories. For each extra constraint
    specified, looks in the fits header to check whether constraint is
    satisfied. If file is good, keep it in temp directory. If file is bad,
    delete it and its .LOG counterpart.
    """
    print(constraints)
    def progress(filename, size, sent):
        sys.stdout.write("%s\'s progress: %.2f%%   \r" % (filename, float(sent) / float(size) * 100))
    rpath = 'R'
    if os.path.exists(rpath):
        while os.path.exists(rpath):
            rpath = rpath + '0'
    os.system('mkdir ' + rpath)
    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.connect('astrometry.mit.edu', username='urop1', password='DEAPS student1')
    with SCPClient(ssh.get_transport(), sanitize=lambda x: x, progress=progress) as scp:
        scp.get('/ast2/data/' + constraints['telescope'] + '/' + constraints['date'][:4] + '/' + constraints['date'].split('/')[0] +
                '/' + constraints['date'].split('/')[0] + '.R.*.gz')
    ssh.close()
    for file in glob.glob(constraints['date'].split('/')[0] + '.R.*.gz'):
        shutil.move(file, rpath + '/')
    filenames = os.listdir(rpath + '/')
    score = 0
    validated = []
    for f in filenames:
        if f[0] == '.':
            continue
        os.system("gunzip '" + rpath + '/' + f + "'")
        for key in constraints.keys():
            if key == 'instrument':
                value = fits.open(rpath + '/' + f[:-3])[0].header['INSTRUME']
                if value != constraints[key]:
                    score = score + 1
            elif key == 'object':
                value = fits.open(rpath + '/' + f[:-3])[0].header['OBJECT']
                if value != constraints[key]:
                    score = score + 1
        if score == 0:
            validated.append(f[-7:-3])
        else:
            print('removing ' + f)
            os.system("rm -rf " + rpath + "/" + f[:-3])
        score = 0
    return validated

def readnumbers():
    """Gets the numbers of the files that satisfied constraints and are sitting in the R/ file"""
    filenames = os.listdir(os.path.join(os.getcwd(),'R'))
    numbers= []
    for name in filenames:
        numbers.append(name[-4:])
    return numbers

def fetch(telescope, date, type, numbers=None): #type is R, Rb, or log
    def progress(filename, size, sent):
        sys.stdout.write("%s\'s progress: %.2f%%   \r" % (filename, float(sent) / float(size) * 100))
    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.connect('astrometry.mit.edu', username='urop1', password='DEAPS student1')
    with SCPClient(ssh.get_transport(), sanitize=lambda x: x, progress=progress) as s:
        if type == 'R':
            if numbers is None:
                s.get('/ast2/data/' + telescope + '/' + date[:4] + '/' + date.split('/')[0] +
                        '/' + date.split('/')[0] + '.R.*.gz')
            else:
                for number in numbers:
                    s.get('/ast2/data/' + telescope + '/' + date[:4] + '/' + date.split('/')[0] +
                            '/' + date.split('/')[0] + '.R.' + number + '.gz')
        elif type == 'Rb':
            if numbers is None:
                s.get('/ast2/data/' + telescope + '/' + date[:4] + '/' + date +
                    '/' + date.split('/')[0] + '.Rb.*')
            else:
                for number in numbers:
                    try:
                        s.get('/ast2/data/' + telescope + '/' + date[:4] + '/' + date +
                                '/' + date.split('/')[0] + '.Rb.'+ number)
                    except scp.SCPException:
                        print("There is no Rb file for " + str(number))
                        os.system('rm -rf R/' + date + '.Rb.' + number)

        elif type == 'log':
            if numbers is None:
                s.get('/ast2/data/' + telescope + '/' + date[:4] + '/' + date +
                    '/' + date.split('/')[0] + '.Rb.*.LOG')
            else:
                for number in numbers:
                    s.get('/ast2/data/' + telescope + '/' + date[:4] + '/' + date +
                            '/' + date.split('/')[0] + '.Rb.' + number + '.LOG')
        else:
            raise Exception('type must be R, Rb, or log')
    ssh.close()
    if type == 'R':
        if os.path.exists(type):
            while os.path.exists(type):
                type = type + '0'
        os.system('mkdir ' + type)
        for file in glob.glob(date.split('/')[0] + '.R.*.gz'):
            shutil.move(file, type +'/')
    elif type == 'Rb':
        if os.path.exists(type):
            while os.path.exists(type):
                type = type + '0'
        os.system('mkdir ' + type)
        for file in glob.glob(date.split('/')[0] + '.Rb.*'):
            shutil.move(file, type + '/')
    elif type == 'log':
        if os.path.exists(type):
            while os.path.exists(type):
                type = type + '0'
        os.system('mkdir ' + type)
        for file in glob.glob(date.split('/')[0] + '.Rb.*.LOG'):
            shutil.move(file, type + '/')
    else:
        raise Exception('type must be R, Rb, or log')
    print(str(len(os.listdir(type + '/'))) + ' ' + type + ' files successfully copied to ' + str(os.path.join(os.getcwd(), type)))


def run_cmd(sshClient, command):
    channel = sshClient.get_transport().open_session()
    channel.get_pty()
    channel.exec_command(command)
    out = channel.makefile().read()
    err = channel.makefile_stderr().read()
    returncode = channel.recv_exit_status()
    channel.close()                       # channel is closed, but not the client
    return out, err, returncode


def readcsv(filepath, parse_nulls=True, fdp = False):
    """Reads an Rb file into a pandas dataframe"""
    data = {}
    keys = []
    fdp_data = {'parameter1':[], 'uncertainty1':[],
                'parameter2':[],'uncertainty2':[],'Fit1':[],'Fit2':[]}
    fdp_keys = ['parameter1', 'uncertainty1','parameter2',
                'uncertainty2', 'Fit1','Fit2']
    fdp_chi2 = []
    csv = open(filepath)
    keyindex = 0
    startindex = 0
    endindex = 1
    substring = ""
    linenumber = 0
    csv.readline()
    lines = csv.readlines()
    startline = len(lines)
    i = 0
    while i < len(lines):
        if lines[i][:9] == '#MATCHTOL' or i > startline:
            if lines[i][:9] == '#MATCHTOL':
                startline = i
                i = i + 2
            if lines[i][:19] == '#PlateFitChiSquared':
                startindex = 14
                endindex = startindex + 1
                keyindex = 0
                while endindex <= len(lines[i]):
                    substring = lines[i][startindex:endindex]
                    if substring[-1] == "." or (substring[0] == "-" and len(substring) == 1):
                        endindex = endindex + 1
                        substring = lines[i][startindex:endindex]
                    try:
                        number = float(substring)
                        if substring[-1] == " " or substring[-1] == "\t" or substring[-1] == "\n" or endindex == len(
                                lines):
                            fdp_chi2.append(number)
                            keyindex = keyindex + 1
                            startindex = endindex
                        endindex = endindex + 1
                    except ValueError:
                        startindex = endindex
                        endindex = endindex + 1
            if len(fdp_data['parameter1']) == 14:
                break
            startindex = 9
            endindex = startindex+1
            keyindex = 0
            while endindex <= len(lines[i]):
                substring = lines[i][startindex:endindex]
                if substring[-1] == "." or (substring[0] == "-" and len(substring) == 1) or substring[-1] == 'e':
                    endindex = endindex + 1
                    if substring[-1] == 'e':
                        endindex = endindex + 1
                    substring = lines[i][startindex:endindex]
                    #if linenumber == 46: #DEBUG
                        #print("line:", linenumber, " start:", startindex, " end:", endindex, ' substring:', substring, ' keyindex:', keyindex) #DEBUG
                try:
                    number = float(substring)
                    if substring[-1] == " " or substring[-1] == "\t" or substring[-1] == "\n" or endindex == len(lines):
                        #if linenumber == 46: #DEBUG
                        #print("appending " + substring + " to " + fdp_keys[keyindex]) #DEBUG
                        fdp_data[fdp_keys[keyindex]].append(number)
                        keyindex = keyindex + 1
                        startindex = endindex
                    endindex = endindex + 1
                except ValueError:
                    startindex = endindex
                    endindex = endindex + 1
        i = i+1
    csv.close()
    ###############################################################################
    csv = open(filepath)
    for line in csv:
        linenumber = linenumber + 1
        if line[:3] == '#1 ':
            startindex = 3
            endindex = 4
            while endindex <= len(line):
                substring = line[startindex:endindex]
                if substring[-1] == " " or substring[-1] == '\t' or endindex == len(line):
                    if substring.count(' ') != len(substring):
                        substring = substring.strip(' ')
                        substring = substring.strip('\n')
                        keys.append(substring)
                        startindex= endindex
                endindex= endindex + 1
            for key in keys:
                data[key] = []

        if line[0] == "#" or line [:2] == '\n':
            continue
        startindex = 0
        endindex = 1
        keyindex = 0
        while endindex <= len(line):
            substring = line[startindex:endindex]
            if substring[-1] == "." or (substring[0] == "-" and len(substring) == 1):
                endindex = endindex + 1
                substring = line[startindex:endindex]
            #if linenumber == 234: #DEBUG
                #print("line:", linenumber, " start:", startindex, " end:", endindex, ' substring:', substring, ' keyindex:', keyindex) #DEBUG
            try:
                number = float(substring)
                if substring[-1] == " " or substring[-1] == "\t" or substring[-1] == "\n" or endindex == len(line):
                    #if linenumber == 234: #DEBUG
                        #print("appending " + substring + " to " + keys[keyindex]) #DEBUG
                    data[keys[keyindex]].append(number)
                    keyindex = keyindex + 1
                    startindex = endindex
                endindex = endindex + 1
            except ValueError:
                startindex = endindex
                endindex = endindex + 1

    if parse_nulls:
        data = pd.DataFrame(data)
        data = data[[int(x) > 0 for x in data['RefID']]]
        data = data.reset_index()
    if not fdp:
        return data
    else:
        return data,(fdp_data, fdp_chi2)

#DEBUG READCSV#
#data, fdp = readcsv('mitastrometry/SARAS/20160611/20160611.Rb.150',fdp=True)
#data = pd.DataFrame(data)
#print(data, fdp)

def fetchlmi(date, object, filetype):
    """
    wrapper function for fetch() and satisfies_constraints() on DCT/LMI data
    scp's the data you want and deletes what you don't want.
    works by taking all the R files from the date dir, unzipping them,
    putting them in an R directory, and using the fits headers to delete
    everything that doesn't satisfy object constraint.
    Passes the numbers that satisfy object constraint to fetch(), which gets only
    the Rb/log files you want and puts them in an Rb/log folder.
    Run this function in the directory you want the R and Rb/log files of your
    target object in.
    It would be more efficient to use grep to find the files that satisfy
    constraints, but for reasons I can't recall I used this method instead
    ARGUMENTS
    ---------
    date (str) : Target date eg. '20190901a'
    object (str): Target object eg. '29P'
    filetype (str):'Rb' or 'log'
    """
    telescope = 'DCT/LMI'
    fetch(telescope, date + '/linearFDP2x2_GDR2', filetype,
        numbers = satisfies_constraints({'telescope': telescope,
                                        'date':date,
                                        'object':object}) )


#eg. fetchlmi('20191021b', '29P', 'Rb')
