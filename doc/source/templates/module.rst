:mod:`{{module}}`.{{objname}}
{{ underline }}==============

..
    .. currentmodule:: {{ module }}

.. automodule:: {{ module }}.{{ objname }}

.. autosummary::
    :toctree: ./
    :template: class.rst

    {% for cls in classes %}
    {{ cls }}
    {% endfor %}


.. raw:: html

    <div class="clearer"></div>
