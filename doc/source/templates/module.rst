:mod:`{{module}}`.{{objname}}
{{ underline }}==============

..
    .. currentmodule:: {{ module }}

.. automodule:: {{ module }}.{{ objname }}


{% if classes %}

Classes
-------

.. autosummary::
    :toctree: ./
    :template: class.rst

    {% for cls in classes %}
    {{ cls }}
    {% endfor %}

{% endif %}

{% if functions %}

Functions
---------

.. autosummary::
    :toctree: ./

    {% for f in functions %}
    {{ f }}
    {% endfor %}

{% endif %}


.. raw:: html

    <div class="clearer"></div>
