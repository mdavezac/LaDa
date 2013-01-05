import sys
sys.path = ["@CMAKE_CURRENT_BINARY_DIR@/..", "@CMAKE_CURRENT_BINARY_DIR@/.."] + sys.path
from pylada import error
import exception_@TYPE@

exception_@TYPE@.nothrow()
try: exception_@TYPE@.dothrow_nomessage()
except error.@TYPE@: pass
try: exception_@TYPE@.dothrow_message()
except error.@TYPE@ as e:  
  assert str(e).find("This is a message.") != -1
try: exception_@TYPE@.dopythrow_message()
except error.@TYPE@ as e: 
  assert str(e).find("This is another message.") != -1

try: raise error.@TYPE@("Whatever")
except error.@TYPE@ as e: 
  assert str(e).find("Whatever") != -1
try: raise error.@TYPE@("Whatever")
except @TYPE@ as e: 
  assert str(e).find("Whatever") != -1
