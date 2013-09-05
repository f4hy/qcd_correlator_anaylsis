#!/usr/bin/env python
import sys
import os
import re

bold = "\033[1;37m"
reset = "\033[0;0m"
queues = ("green", "blue", "cyan", "magenta", "red")


def choices(c):
    return dict(enumerate(c))


def selectchoices(c, default=None, startnum=0):

    if default and default not in c:
        raise Exception("default choice %s not in choices" % str(default))

    if default:
        prompt = "please select\n"
        for x in enumerate(c, start=startnum):
            if x[1] == default:
                prompt += bold
                prompt += "[%d] %s default\n" % x
                prompt += reset
            else:
                prompt += "[%d] %s \n" % x
        # prompt = ("please select\n" + "\n".join(["[%d] %s default" % x if x[1] == default
        #                                         else "[%d] %s " % x
        #                                         for x in enumerate(c)]) + "\n")
    else:
        prompt = ("please select\n" + "\n".join(["[%d] %s" % x for x in enumerate(c)]) + "\n:")

    dchoices = dict(enumerate(c, start=startnum))
    while True:
        sys.stdout.write(prompt)
        choice = raw_input().lower()
        if choice == "exit":
            raise Exception('exit')
        elif default is not None and choice == '':
            return default
        elif choice in c:
            return choice
        elif choice in dchoices:
            return dchoices[choice]
        else:
            try:
                intchoice = int(choice)
            except ValueError:
                sys.stdout.write("Invalid selection %s \n" % choice)
            else:
                if intchoice in dchoices:
                    return dchoices[intchoice]
                else:
                    sys.stdout.write("Invalid: selection %d out of range \n" % intchoice)
    # End select choices


def askqueue(defaultqueue="red"):
    queue = selectchoices(queues, default=defaultqueue)
    sys.stdout.write("selected %s queue\n" % queue)
    return queue
    #return selectchoices(queues)


def askrange(start=1, end=1, defaultval=1):
    defaultval = min(max(start, defaultval), end)     # make sure default is in the range
    return selectchoices(list(range(start, end + 1)), default=defaultval, startnum=start)


def askdir(description, default=os.getcwd()):
    prompt = "Enter directory for %s \n (%s):" % (description, default)

    while True:
        sys.stdout.write(prompt)
        userdir = raw_input()
        if userdir == "exit":
            raise Exception('exit')
        elif userdir == '' and os.path.isdir(default):
            return default
        elif os.path.isdir(default):
            return default
        elif os.path.isdir(userdir):
            return userdir
        else:
            sys.stdout.write("INVALID direcotry %s\n" % userdir)


def askoperator(description, default=""):
    prompt = "Enter operator for %s \n (%s):" % (description, "operator, exit, or skip")
    while True:
        sys.stdout.write(prompt)
        userstring = raw_input()
        operator_pattern = re.compile("[A-z]+.*\_[0-9]+")
        if userstring == "exit":
            raise Exception('exit')
        elif userstring == "skip":
            sys.stdout.write("Skipping")
            return None
        elif operator_pattern.match(userstring):
            return userstring
        else:
            sys.stdout.write("INVALID entry %s\n" % userstring)


def askstring(description, default=""):
    prompt = "Enter directory for %s \n (%s):" % (description, default)
    while True:
        sys.stdout.write(prompt)
        userstring = raw_input()
        if userstring == "exit":
            raise Exception('exit')
        elif userstring == '' and default:
            return default
        elif userstring:
            return userstring
        else:
            sys.stdout.write("INVALID entry %s\n" % userstring)


def askyesno(description, default=True):
    if default:
        prompt = "%s [%sY%s]/n " % (description, bold, reset)
    else:
        prompt = "%s y/[%sN%s] " % (description, bold, reset)
    while True:
        ans = raw_input(prompt)
        if not ans:
            return default
        if ans.lower() not in ['y', 'yes', 'n', 'no']:
            print 'please enter y or n.'
            continue
        if ans.lower() == 'y' or ans == 'yes':
            return True
        if ans.lower() == 'n' or ans == 'no':
            return False


def readgeom(default=(1, 1, 1, 1)):
    while True:
        ans = raw_input("Input gemoetry ('%d %d %d %d')" % tuple(default))
        if not ans:
            return default
        try:
            x, y, z, t = ans.split(' ')
            print x, y, z, t
            print ans.split(' ')
            if not all(s.isdigit() for s in ans.split(' ')):
                print "not valid numbers"
                continue
        except ValueError:
            print 'invalid form for geom must be like "1 1 1 2"'
            continue
        return [int(i) for i in ans.split(' ')]


if __name__ == '__main__':
    print askqueue()
    # selected = selectchoices(queues)
    # print "selected thing was %s" % selected
