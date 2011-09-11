import sys
sys.path = ["@CMAKE_CURRENT_BINARY_DIR@/..", "@CMAKE_CURRENT_BINARY_DIR@/.."] + sys.path
import error
import exception_@TYPE@

exception_@TYPE@.nothrow()
try: exception_@TYPE@.dothrow_nomessage()
except error.@TYPE@: pass
try: exception_@TYPE@.dothrow_message()
except error.@TYPE@ as e: 
  print str(e)
  assert str(e).find("This is a message.") != -1
