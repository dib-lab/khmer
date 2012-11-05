.. labibi documentation master file, created by
   sphinx-quickstart on Sun Nov  4 10:10:29 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the labibi demo site!
================================

.. toctree::
   :maxdepth: 2

Labibi is a base package to use for documentation and Web sites for
other projects of mine.

For now, I've made it fairly easy to post sites to `ReadTheDocs
<http://readthedocs.org>`__ that enable Google Analytics, Disqus commenting,
and easy source file editing.  You can check it out at the `labibi demo
site <http://labibi.readthedocs.org/en/latest/>`__.

This is directly based off of Mikko Ohtamaa's excellent work on `the
Plone documentation
<http://opensourcehacker.com/2012/01/08/readthedocs-org-github-edit-backlink-and-short-history-of-plone-documentation/>`__.

A brief HOWTO do this for your own ReadTheDocs site:

  0. `Get started with ReadTheDocs <https://docs.readthedocs.org/en/latest/getting_started.html>`__.

  1. Create a _static/ directory and put `labibi.css <https://raw.github.com/ctb/labibi/master/_static/labibi.css>`__ and `labibi.js <https://raw.github.com/ctb/labibi/master/_static/labibi.js>`__ in it.

  2. Put "html_style = 'labibi.css'" in your conf.py

  3. Create a _templates/ directory and put `page.html <https://raw.github.com/ctb/labibi/master/_templates/page.html>`__ in there.

  4. Edit 'page.html' to set your google analytics, disqus, and github info.

For now, you can't disable the editing functionality, but if you
delete the google_analytics and disqus_shortname it should disable
that functionality on your site.

The labibi source code is `here <https://github.com/ctb/labibi>`__ for
your checking-out pleasure.

I've put up a short blog post `here <http://ivory.idyll.org/blog/rtd-comments-and-editing.html>`__.

--titus


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

