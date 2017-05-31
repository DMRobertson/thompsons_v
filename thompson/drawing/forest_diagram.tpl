{% if standalone %}
\documentclass[tikz]{standalone}
\usetikzlibrary{graphs, graphdrawing, positioning, calc, backgrounds}
\usegdlibrary{trees}
\begin{document}
{% endif %}
{% if include_styles: %}
%from http://tex.stackexchange.com/a/47004/82389
\pgfdeclarelayer{back}
\pgfsetlayers{back,main}
\pgfkeys{ 
  /tikz/on layer/.code={
    \pgfonlayer{{'{#1}'}}\begingroup
    \aftergroup\endpgfonlayer
    \aftergroup\endgroup
  }
}

\tikzset{
	domain tree/.style = { name = domain },
	range tree horiz LTR/.style = {
		name = range,
		right=2.4em of domain.north east,
		anchor=north west,
	},
	range tree horiz RTL/.style = {
		name = range,
		left=2.4em of domain.north west,
		anchor=north east,
	},
	range tree vert/.style = {
		name = range,
		below=5em of domain.south,
		anchor=north,
	},
	connecting arrow/.style = {
		->, semithick
	},
	component sep=3em,
	component/.style   = {dotted},
	spine/.style       = {
		component,
		preaction={
			solid,
			draw=red!15,
			line width=1mm,
			line cap=round,
			on layer=back
		}
	},
	repatt/.style      = {red},
	repatt label/.style= {red, fill=white},
	every label/.style = {fill=white, inner sep=0pt, outer sep=3.5pt}
}

\tikzgraphsset{
	every graph/.style = {
		extended binary tree layout,
		simple,
		empty nodes,
		edges = { line cap = round },
		nodes = {
			inner sep = 0, 
			minimum size = 0,
		},
		level distance = 2.25em,
		level sep = 0,
		sibling distance = 2.25em,
		sibling sep = 0,
		significant sep = 0.1em,		
	}
}
{% endif %}
\begin{tikzpicture}
{% macro basis_template(basis, style) -%}
	\node [{{style}}] {
		\tikz\graph{
			{% for word, label, highlight in basis %}
				{#{word}} {{label}} {{highlight}#}
			{{ write_word(word, label, highlight) }}
		{% endfor %}
		}; %end \tikz\graph
	}; %end \node
{%- endmacro %}
	{{ basis_template(domain, 'domain tree') }}
{% if horiz %}
	{% if LTR %}
	{{ basis_template(range, 'range tree horiz LTR') }}
	{% else %}
	{{ basis_template(range, 'range tree horiz RTL') }}
	{% endif %}
{% else %}
 	{{ basis_template(range, 'range tree vert' ) }}
{% endif %}
\draw[connecting arrow]
	{% if horiz and LTR %}
	let \p1=(domain.east), \p2=(range.west), \n1={max(\y1,\y2)} in
		(\x1, \n1) -- node[above] { {{name}} } (\x2, \n1);
	{% elif horiz %}
	let \p1=(domain.west), \p2=(range.east), \n1={max(\y1,\y2)} in
		(\x1, \n1) -- node[above] { {{name}} } (\x2, \n1);
	{% else %}
		(domain) -- node[left] { {{name}} } (range);
	{% endif %}
\end{tikzpicture}
{% if standalone %}
\end{document}
{% endif %}