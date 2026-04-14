      program helloworld
      implicit none
      character*13 hello_string
      hello_string = "Hello, world!" !THIS IS ONE WAY TO PRINT A STRING
      write (*,*) hello_string
      print *, "Hello, World!" !THIS IS ANOTHER, THOUGH NOT AS MODULAR
      end program helloworld