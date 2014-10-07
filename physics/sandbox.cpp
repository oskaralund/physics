#include "window.h"
#include "rod_sandbox.h"
#include "pbd_sandbox.h"

#include <iostream>
#include <memory>

int main()
{
  std::cout << "1. Rod sandbox" << std::endl;
  std::cout << "2. PBD sandbox" << std::endl;

  std::unique_ptr<Window> window = nullptr;

  switch (getchar())
  {
  case '1':
    window = std::make_unique<RodSandbox>();
    break;
  case '2':
    window = std::make_unique<PBDSandbox>();
    break;
  }

  if (window)
    window->OpenWindow();

  return 0;
}
