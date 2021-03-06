-- sketch-domain.lua

s = Sketch:new{renderer="svg", projection="xyortho", canvas_mm={0.0,0.0,120.0,120.0}}
s:set{viewport={-0.3,-0.2,1.2,1.2}}

s:start{file_name="heatshield.svg"}
s:text{point=Vector3:new{x=0.5,y=1.08},
       text="flow domain for heatshield",
       font_size=16}

-- for flow domain
s:set{line_width=0.1, fill_colour="green"}
s:render{surf=quad0}; s:render{surf=quad1}

s:set{line_width=0.5}; s:render{path=bc} 
s:dotlabel{point=a, label="a"}; s:dotlabel{point=b, label="b"}
s:dotlabel{point=c, label="c"}; s:dotlabel{point=d, label="d"}
s:dotlabel{point=e, label="e"}; s:dotlabel{point=f, label="f"}
s:dotlabel{point=g, label="g"}; s:dotlabel{point=h, label="h"}
s:dotlabel{point=i, label="i"}; s:dotlabel{point=j, label="j"}
s:dotlabel{point=k, label="k"}; s:dotlabel{point=l, label="l"}


s:set{line_width=0.3} -- for drawing rules
s:rule{direction="x", vmin=0.0, vmax=1.0, vtic=0.2,
       anchor_point=Vector3:new{x=0,y=-0.05},
       tic_mark_size=0.02, number_format="%.1f",
       text_offset=0.06, font_size=10}
s:text{point=Vector3:new{x=0.5,y=-0.15}, text="x", font_size=12}
s:rule{direction="y", vmin=0.0, vmax=1.0, vtic=0.2,
       anchor_point=Vector3:new{x=-0.05,y=0},
       tic_mark_size=0.02, number_format="%.1f",
       text_offset=0.02, font_size=10}
s:text{point=Vector3:new{x=-0.2,y=0.5}, text="y", font_size=12}

s:finish{}
