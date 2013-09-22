==============
Using 'screen'
==============

:Author: Rosangela Canino-Koning
:Date: June 9, 2011
:Last Updated: July 24, 2013

Persistent Sessions
-------------------

Screen is a window manager for terminal sessions. Screen allows you to
run a terminal session, and then disconnect from the computer, and be
able to return to the session at a later date.

To start screen, you run the screen command with a few options::

  screen -S <sessionname>
 
Where *sessionname* is any meaningful or descriptive title for your screen 
session. This creates an independent terminal session, and connects you to it. 

Most commands within screen are composed of a prefix key-stroke,
followed by a command character. By default, the prefix is Ctrl-A. In
this tutorial Ctrl-A will represented by "C-a".

Let's try a few screen commands.

To disconnect from the session (while leaving it running!), type::

  C-a d

This session will remain active until you choose to end it, or you
reboot the computer. You can at this point safely disconnect from SSH,
and the screen session will continue to run.

To reconnect to the session, make sure you're logged into the UNIX machine,
and type::

  screen -r

To illustrate managing multiple screen session, disconnect from the current 
session, and create a new session with a second name.::

  C-a d
  screen -S <secondsessionname>
 
Disconnect from the second session, and then list the available sessions::

  C-a d
  screen -list

Note, typing *screen -r* with multiple active screen sessions will display
the same information.

To reconnect to the first session, include its name after the -r.::

  screen -r <sessionname>
