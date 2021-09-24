
.. _contribute:

Contributing
============

Thanks for contributing! If you would like to suggest large-scale refactors, feel free submit a `GitHub issue <https://github.com/RomeroLab/schemarecomb/issues>`_. Otherwise, you can follow this guide for writing and submitting features for review. (There's a few equivalent ways to do each of the steps.)

We'll use the Fork-and-Branch Workflow. Here's more information on it: `https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/`.


Requirements: schemarecomb requires Python 3.9. This guide assumes that the 'python' command calls Python 3.9.


1. `Fork the repository <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_. 

2. Install `Poetry <https://python-poetry.org/>`_.

3. In your workspace directory (for me this is "~/Programming"), make a new virtual environment, cd into it, and activate it. Poetry will then install everything to the virtual environment, rather than globally.

.. code:: bash

   cd <workspace directory>
   python -m venv schemarecomb_env
   cd schemarecomb_env
   . bin/activate

4. Clone your forked schemarecomb repo and install with Poetry. This will install the development requirements to the virtual environment, so as long as you have schemarecomb_env activated, you should have access to the following commands.

.. code:: bash

   git clone https://github.com/<GitHub username>/schemarecomb
   cd schemarecomb
   poetry install

5. Make a new branch. You should make the branch name something informative, like "feature/local-muscle" if you were adding the ability to run MUSCLE locally.

.. code:: bash

   git branch -b <branch name>

6. Add the desired source code changes and tests.

7. Run the pytest, mypy, and flake8. If there are print-outs for any of these, fix the offending code and rerun. Note these commands are ran from the schemarecomb_env/schemarecomb directory.

.. code:: bash

   pytest
   flake8 src/schemarecomb
   mypy src/schemarecomb
   mypy tests/unit

8. Make a new commit with the source code and test files you modified. If you're unfamiliar with git, `here <https://git-scm.com/book/en/v2/Git-Basics-Recording-Changes-to-the-Repository>`_ is a good place to start learning.

9. Push the changes to your forked GitHub repo.

.. code:: bash

   git push origin <branch name>

10. `Open a pull request. <https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork>`_. A package maintainer will review the pull request and ask for changes if necessary. Otherwise, they'll accept your contribution.
