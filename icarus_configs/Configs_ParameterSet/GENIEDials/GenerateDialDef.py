import os
import knoblist

# NSigmas

knobNames = knoblist.GENIEMultisigmaKnobNames

for name in knobNames:

  IsEnabled = 'true'

  output = '''- parameterName: "%s"
  isEnabled: %s
  dialSetDefinitions:
    - dialType: Spline
      minimumSplineResponse: 0
      dialLeafName: "%s"
      applyCondition: "[IsData]==0"
'''%(name, IsEnabled, name)
  print(output)

# Morph

knobNames = knoblist.GENIEMorphKnobNames

UseMirror = 'true'
for name in knobNames:
  output = '''- parameterName: "%s"
  isEnabled: %s
  dialSetDefinitions:
    - dialType: Spline
      minimumSplineResponse: 0
      dialLeafName: "%s"
      useMirrorDial: true
      mirrorLowEdge: -1
      mirrorHighEdge: 1
      applyCondition: "[IsData]==0"
'''%(name, UseMirror, name)
  print(output)

