testsuite.nonnorm <- defineTestSuite("nonnorm",
                     dir=file.path("/Users/testsuite"),testFileRegexp = "^runit.+\\.[rR]$",
                     testFuncRegexp = "^test.+",
                     rngKind = "Marsaglia-Multicarry",
                     rngNormalKind = "Kinderman-Ramage")
testResult <- runTestSuite(testsuite.nonnorm)
printTextProtocol(testResult)
