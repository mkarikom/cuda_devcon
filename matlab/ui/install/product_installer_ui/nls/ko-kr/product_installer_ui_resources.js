/* Copyright 2020-2021 The MathWorks, Inc. */
define({
    
        dialogTitle: 'MATLAB 인스톨러',
        installButtonLabel: '설치 시작',
        downloadButtonLabel: '다운로드 시작',
        forwardButtonLabel: '다음',
        cancelButtonLabel: '취소',
        finishButtonLabel: '닫기',
        guiltForwardButtonLabel: '학생입니다.',
        closeLabel: '확인',
        continueLabel: '계속',
        ignoreLabel: '무시',
        licenseAgreementTitle: 'MathWorks 라이선스 계약',
        patentsAndTrademarksTitle: '저작권, 상표 및 특허',
        fikWorkflowLabel: '파일 설치 키를 사용하여 설치',
        proxyTitle: '프록시 인증',
        proxyUser: '프록시 서버 사용자 ID:',
        proxyPassword: '프록시 서버 비밀번호:',
        downloadOnlyWorkflowLabel: '설치하지 않고 다운로드',
        licenseManagerWorkflowLabel: '네트워크 라이선스 매니저 설치',
        fikLabel: '파일 설치 키 입력',
        fikContextHelp: '파일 설치 키는 사용자가 설치할 수 있는 제품을 식별하고, 인터넷에 연결하지 않고도 제품을 설치할 수 있게 해 줍니다.  이 방법으로 설치하려면 라이선스 파일도 있어야 합니다. MathWorks 라이선스 센터나 라이선스 관리자로부터 파일 설치 키와 라이선스 파일을 모두 받을 수 있습니다.',
        downloadFolderContextHelp: '나중에 한 대 이상의 컴퓨터에 제품을 설치할 수 있도록 제품을 다운로드하지만 지금 설치하지는 않습니다.',
        activationKeyContextHelp: '활성화 키는 사용자를 다른 라이선스에 연결합니다. 라이선스 관리자로부터 활성화 키를 받을 수 있습니다.',
        proxyContextHelp: '인터넷에 액세스할 때 프록시 서버를 사용하도록 컴퓨터가 구성되었습니다. 표준 온라인 설치를 계속하려면 사용자 이름과 비밀번호를 입력하십시오. 프록시 서버 자격 증명에 대해서는 시스템 관리자에게 문의하십시오.',

        onlineWorkflowButton: '표준 설치를 수행합니다.',
        downloadOnlyWorkflowButton: '설치하지 않고 다운로드합니다.',
        fikWorkflowButton: '파일 설치 키가 있습니다.',
        licenseManagerWorkflowButton: '네트워크 라이선스 매니저를 설치합니다.',

        AKToolTipHelpWithAssociation: '활성화 키는 사용자를 다른 라이선스에 연결합니다. 라이선스 관리자로부터 활성화 키를 받을 수 있습니다.',
        AKToolTipHelpNoAssociation: '활성화 키는 사용자를 라이선스에 연결합니다. 라이선스 관리자로부터 활성화 키를 받을 수 있습니다.',
        downloadOnlyToolTipHelp: '나중에 한 대 이상의 컴퓨터에 제품을 설치할 수 있도록 제품을 다운로드하지만 지금 설치하지는 않습니다.',
        selectLicenseFileForLMToolTipHelp: 'MathWorks 라이선스 센터에서 서버용 라이선스 파일을 생성합니다.',
        selectLicenseFileForToolTipHelp: 'MathWorks 라이선스 센터 또는 라이선스 관리자로부터 라이선스 파일을 받습니다.',
        licenseManagerToolTipHelp: '이 워크플로는 라이선스 관리자만 이용 가능합니다. 동시사용자(Concurrent) 라이선스와 명명된 네트워크 사용자(Network Named User) 라이선스를 사용하려면 서버에 네트워크 라이선스 매니저를 설치해야 합니다.',
        akLabel: '활성화 키 입력',
        installationFolderLabel: '대상 폴더 선택',
        downloadFolderLabel: '대상 폴더 선택',
        platformSelectionLabel: '대상 플랫폼 선택',
        browseButtonLabel: '찾아보기',
        restoreDefault: '디폴트 값 복원',
        productTableLabel: '제품 선택',
        productTableLabelMinimal: '제품 선택(권장되는 제품이 미리 선택되어 있음)',
        productColumnHeader: '모두 선택',
        advancedOptionText: '고급 옵션',
        optionsTitle: '옵션 선택',
        shortcutLabel: '바탕 화면에 바로 가기 추가',
        symbolicLinkLabel: '다음에 MATLAB 스크립트의 기호화된 링크 만들기:',
        userExperienceLabel: '더 개선된 MATLAB이 나올 수 있도록 사용자 경험 정보를 MathWorks로 전송',
        finalLabel: '설치 완료',
        finalLabelDownloadOnly: '다운로드 완료',
        fikNavTitle: '파일 설치 키',
        folderSelectionNavTitle: '설치 폴더',
        productSelectionNavTitle: '제품',
        optionsNavTitle: '옵션',
        confirmationNavTitle: '확인',
        licenseSelectionTitle: '라이선스 파일 선택',
        licenseSelectionFolderLabel: '(라이선스 파일의 전체 경로를 파일 이름을 포함하여 입력)',
        entitlementPanelTitle: '라이선스 선택',
        activationPanelTitle: '활성화 키 입력',
        activationPanelTitleA2AK: '라이선스 활성화 키 입력',
        A2AKPanelContextHelp: '이 라이선스에 액세스하려면 활성화 키가 필요합니다. 라이선스 관리자로부터 활성화 키를 받을 수 있습니다.',
        entitlementTableTitle: '라이선스:',
        entitlementTableLicenseCol: '라이선스',
        entitlementTableLabelCol: '레이블',
        entitlementTableOptionCol: '라이선스 사용 및 옵션',
        activationKeyTitle: '활성화 키 입력:',
        activationModePanelTitle: '인증 옵션 선택',
        activationModeRadioLabel: '지금 이 컴퓨터 인증',
        activationModeRadioText: '이 컴퓨터에서 소프트웨어를 활성화합니다. 인터넷 연결 없이 소프트웨어를 실행할 수 있습니다.',
        activationModeAuthRadioLabel: '소프트웨어가 시작될 때마다 인증(로그인)',
        activationModeAuthRadioText: '소프트웨어를 실행하려면 인터넷 연결이 필요합니다.',
        configureService: '서비스로 구성',
        configureServiceText: 'Windows 시스템에서 인스톨러는 시스템 시작 시 라이선스 매니저가 자동으로 시작되도록 구성합니다. 이 디폴트 구성을 수락한 경우 라이선스 매니저를 시작하는 가장 쉬운 방법은 라이선스 매니저를 설치한 컴퓨터를 다시 시작하는 것입니다.',
        confirmUserTitle: '사용자 확인',
        confirmUserNameLabel: '이름',
        confirmUserEmailLabel: '이메일',
        confirmWinUserNameLabel: 'Windows 사용자 이름',
        editWinUserNameLabel: '편집',
        ConfirmLinuxUserName: 'glnxa64 로그인 이름',
        ConfirmMacUserName: 'maci 로그인 이름',
        linuxOrMacUserNameLabel: '로그인 이름',
        firstName: '이름',
        lastName: '성',
        getUserInfoBtnLabel: '본인이 이 소프트웨어를 사용함',
        openDownloadFolder: '다운로드 폴더 열기',
        activateMATLAB: '소프트웨어 활성화',
        activateNote: '참고: 소프트웨어를 사용하려면 먼저 활성화해야 합니다.',
        provideUserInfoTitle: '사용자 정보 제공',
        provideUserInfoText: '한 사람만 이 라이선스를 사용할 수 있습니다. 라이선스가 부여된 최종 사용자를 지정하십시오.',
        provideUserInfo: '(필요한 경우 이 사용자의 MathWorks 계정이 만들어짐)',
        productDependecyDialogTitle: '제품 종속성',
        productDependecyDialogAddBtnLabel: '추가',
        productDependecyDialogDoNotAddBtnLabel: "추가 안 함",
        productDependecyDialogCancelBtnLabel: '취소',
        student: '예, 학생입니다.',
        notStudent: '아니요, 학생이 아닙니다.',
        studentConfirmTitle: '학생 이용 정책',
        studentConfirmMessage: '이것은 학생 라이선스입니다. 귀하는 학위를 수여하는 고등 교육 기관에 재학 중인 학생이거나 고등학교/예비 대학의 학생 또는 교사여야 합니다. 이 라이선스는 상업적인 용도나 다른 용도로는 사용할 수 없습니다. MathWorks는 학생을 위한 특별 서비스로 이 라이선스를 제공하는 것이며, 라이선스가 남용되지 않도록 주의해 줄 것을 부탁드립니다.\n\n귀하는 학습 과제에 이 소프트웨어를 사용하는 학교, 대학 또는 대학교의 학생입니까?\n\n',
        yes: '예',
        no: '아니요',
        userNameModifyTitle: 'Windows 사용자 이름 경고',
        userNameModifyMessage1: '소프트웨어가 다음과 다른 사용자 이름으로 실행될 경우에만 변경: ${0}',
        userNameModifyMessage2: ' \n \n사용자 이름에 대한 자세한 내용은 <a href="https://www.mathworks.com/pi_unc_mpi_${0}_${1}" target=\"_blank\">MATLAB Answer</a>를 참조하십시오.',
        loginNameModifyTitle: '로그인 이름 경고',
        loginNameModifyMessage1: '소프트웨어가 다음과 다른 로그인 이름으로 실행될 경우에만 변경: ${0}',
        loginNameModifyMessage2: '\n \n로그인 이름에 대한 자세한 내용은 <a href="https://www.mathworks.com/pi_unc_mpi_${0}_${1}" target=\"_blank\">MATLAB Answer</a>를 참조하십시오',
        userNameRequired: '사용자 이름이 필요함',
        confirmationDialogTitle: '인스톨러 중지',
        confirmationDialogText: '지금 중지하면 설치를 처음부터 다시 시작해야 합니다.<br><br>지금 중지하시겠습니까?',
        mySelectionButton: '내가 선택한 항목 사용',
        recommendedSelectionButton: '권장 항목 사용',
        radioButtonDescriptionForLicenseAgreement: '이 라이선스 계약 조건에 동의하십니까?',
        radioButtonsYesText: '예',
        radioButtonsNoText: '아니요',
        signInButtonText: '로그인',
        connectionErrorTitle: '연결 오류',
        connectionErrorMsg: '응용 프로그램을 MathWorks에 연결할 수 없습니다.<br>이 문제를 해결하려면 <a href="http://www.mathworks.com/pi_dlerr_mpi " target="_blank">고객 지원팀</a>에 문의하십시오.',
        confirmationPanelTitle: '선택 사항 확인',
        confirmDestinationFolderLabel: '대상 폴더',
        confirmPlatformsLabel: '플랫폼',
        confirmProductsLabel: '제품',
        confirmLicensingLabel: '라이선싱',
        installationNotesTitle: '<b>설치하려면 추가적인 구성 단계가 필요할 수 있습니다.</b><br><br>',
        licensingNavTitle: '라이선싱',
        destinationNavTitle: '대상 폴더',
        platformsNavTitle: '플랫폼',
        genericErrorTitle: '예기치 않은 문제가 발생했습니다.',
        genericErrorMessage: '이 문제를 해결하려면 <a href="https://www.mathworks.com/support" target="_blank">기술 지원팀</a>에 문의하십시오.',
        closeButtonLabel: '닫기',
        mathworksTrademarksText: 'MATLAB and Simulink are registered trademarks of The MathWorks, Inc. Please see mathworks.com/trademarks for a list of additional trademarks. Other product or brand names may be trademarks or registered trademarks of their respective holders.',
        mathworksPatentsText: 'MathWorks products are protected by patents (see mathworks.com/patents) and copyright laws. Any unauthorized use, reproduction, or distribution may result in civil and criminal penalties.',
        symbolicLinkFolderLabel: '기호화된 링크 폴더 선택',
        noArchiveErrorMsg: '설치할 제품 파일을 찾을 수 없습니다. 도움이 필요하면 <a href="http://www.mathworks.com/pi_miferr_mpi " target="_blank">MATLAB Answer</a>를 참조하십시오.',
        noArchiveErrorTitle: '설치 파일이 누락됨',
        noInternetNoArchiveErrorMsg1: '문제는 다음 중 하나일 수 있습니다:',
        noInternetNoArchiveErrorMsg2: '<li>인스톨러가 MathWorks에 연결할 수 없습니다. 이는 방화벽 또는 바이러스 백신 소프트웨어로 인한 것일 수 있습니다.  자세한 내용은 <a href="http://www.mathworks.com/pi_icerr_mpi " target="_blank">MATLAB Answer</a>를 참조하십시오.</li>',
        noInternetNoArchiveErrorMsg3: '<li>오프라인 설치에 필요한 제품 파일이 누락되었습니다. 자세한 내용은 <a href="http://www.mathworks.com/pi_miferr_mpi " target="_blank">MATLAB Answer</a>를 참조하십시오.</li>',
        noInternetNoArchiveErrorTitle: '계속할 수 없음'
    });
