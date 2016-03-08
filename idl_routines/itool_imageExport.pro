; this can be thrown into the for loop that creates the iplot
; in this for loop, ihr was incremented
; the naming was done so that the frames would be labled iplot10.jpg, iplot11.jpg, etc.
; iplot export can also do .png, .bmp, and .j2000

  idTool=itgetcurrent(tool=oTool)
  ; Export the window to a JPEG file.
  ; First disable the Export wizard dialog.
  void = oTool->DoSetProperty('Operations/File/Export', 'SHOW_EXECUTION_UI', 0)
  ; Export the entire window.
  void = oTool->DoSetProperty('Operations/File/Export', 'SOURCE', 1)
  void = oTool->DoSetProperty('Operations/File/Export', 'FILENAME', $
                                'iplot'+string(ihr+10,format='(I2)')+'.jpg')
  void = oTool->DoAction('Operations/File/Export')
