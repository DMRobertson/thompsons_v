digraph {
	subgraph periodic {
		{% for orbit in finite_orbits %}
			{% for word in orbit %}
				{{ name(word) }} -> 
			{% endfor %}
			{{ name(orbit[0]) }}
		{% endfor %}
	}
	subgraph AttractorRoots {
		{% for x in attractor_roots %}
		{{ name(x) }}
		{% endfor %}
	}
	subgraph RepellerRoots {
		{% for x in repeller_roots %}
		{{ name(x) }}
		{% endfor %}
	}
	#Flow lines
	{% for chain in chains %}
		{{ name(chain[0]) }}
		{% for word in chain[1:] %}
			{{ ' -> ' + name(word) }}
		{% endfor %}
	{% endfor %}
}
