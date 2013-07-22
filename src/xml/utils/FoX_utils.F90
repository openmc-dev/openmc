module FoX_utils
  
  use fox_m_utils_uuid
  use fox_m_utils_uri

  implicit none
  private

  public :: generate_uuid

  public :: URI
  public :: parseURI
  public :: rebaseURI
  public :: copyURI
  public :: destroyURI
  public :: expressURI
  public :: hasFragment
  public :: hasScheme
  public :: getScheme
  public :: getPath

end module FoX_utils
