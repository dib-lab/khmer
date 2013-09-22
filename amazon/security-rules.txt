========================
Adjusting security rules
========================

Before continuing, you'll need to adjust your security rules so that you
can access your new instance properly.  To do that, go over to Security
Groups on the dashboard, select 'default', and then adjust your security
rules to enable ports 22, 80, and 443 (SSH, HTTP, and HTTPS).

.. image:: images/ec2-security.png
   :width: 90%

Make sure you "Apply rule changes" afterwards.

Then, go to :doc:`log-in-with-ssh-win` or :doc:`log-in-with-ssh-mac`
